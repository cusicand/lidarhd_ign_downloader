#!/usr/bin/env python
# coding=utf-8

import os
import subprocess
import time
import argparse

from pathlib import Path
import wget
import requests
import json
import pdal
from tqdm import tqdm

import geopandas as gpd
import psutil
import concurrent.futures

# percentage of the CPU workload to be used for processing. 
# Warning: keep some CPUs for the OS and other processes(at least 4 CPUs).
CPU_WORKLOAD = 0.6


def process_tile(laz_path, output_path, resolution):
    # This wraps the PDAL processing
    json_pipeline = pdal_json_pipeline(
        str(laz_path), str(output_path), out_type="mean", tr=resolution
    )
    pdal_json_str = json.dumps(json_pipeline)
    pipeline = pdal.Pipeline(pdal_json_str)
    pipeline.execute()


def download_file(url, output_path):
    response = requests.get(url, timeout=10)
    if "content-disposition" in response.headers:
        content_disposition = response.headers["content-disposition"]
        filename = content_disposition.split("filename=")[1]
    else:
        filename = url.split("/")[-1]
    if not os.path.exists(os.path.join(output_path, filename)):
        with open(os.path.join(output_path, filename), mode="wb") as file:
            file.write(response.content)
        print(f"{filename} Downloaded.")
    else:
        print(f"{filename} already exists.")


def pdal_json_pipeline(input_laz_fn, out_tif_fn, out_type="mean", tr=1):
    pdal_json_pipeline = {
        "pipeline": [
            # To do: check if this is a text or posix object
            input_laz_fn,
            {"type": "filters.range", "limits": "Classification[1:2]"},
            {
                "filename": out_tif_fn,
                "gdaldriver": "GTiff",
                "output_type": out_type,
                "resolution": tr,
                "type": "writers.gdal",
            },
        ]
    }
    return pdal_json_pipeline


def getparser():
    args_desc = """
    --------------------------------------------------
    lidarhd downloader --> 'lidarhd_downloader.py' module

    Script to download lidar-hd classified tiles from IGN.
    --------------------------------------------------"""

    # Create an argument parser to get the CLI user arguments
    parser = argparse.ArgumentParser(
        prog="lidarhd_downloader.py.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=args_desc,
    )
    # MANDATORY - Inputs
    parser.add_argument("aoi", type=str, help="path to image aoi.shp file")
    # OPTIONNALS
    parser.add_argument(
        "-out_data",
        "--out_data_path",
        type=str,
        default=None,
        help="Out data-path directory. If not specified, data will be stored in lidar_downloader.py base-path by default.",
    )
    parser.add_argument(
        "-tr",
        "--dem_resolution",
        type=float,
        default=1.0,
        help="resolution of output DEM. Default value : 1.0",
    )
    parser.add_argument(
        "-dtype",
        "--file_data_type",
        type=str,
        default="gtiff",
        choices=["gtiff", "vrt"],
        help="Out data format between GeoTIff and Virtual Dataset (VRT). Default value : GTiff",
    )
    return parser


def main():
    # Get arguments
    parser = getparser()
    args = parser.parse_args()
    if args.out_data_path == None:
        workdir = Path.cwd()
        # print(f"As not out_path have been specified, data will be stored in $Home/lidarhd_ign_downloader/raw_laz_data\n-----")
    else:
        workdir = Path(args.out_data_path)
        if not workdir.exists():
            workdir.mkdir()

        print(f"Data will be stored in {workdir}/raw_laz_data\n-----")

    print(f"Working on: {workdir}")
    extraction_path = workdir.joinpath("raw_laz_data")
    if not extraction_path.exists():
        # Creating directory if not exist
        extraction_path.mkdir()

    tiles_fn = workdir.joinpath("ign_resources", "TA_diff_pkk_lidarhd_classe.shp")

    if not workdir.joinpath("ign_resources").exists():
        workdir.joinpath("ign_resources").mkdir()

    if not tiles_fn.exists():
        print("@INFOs-----\nDownloading tiles shapefile from IGN server")
        wget.download(
            url="https://diffusion-lidarhd-classe.ign.fr/download/lidar/shp/classe",
            out=str(workdir.joinpath("ign_resources")),
        )
        os.system(
            f"unzip {str(workdir.joinpath('ign_resources', 'grille.zip'))} -d {str(workdir.joinpath('ign_resources'))}"
        )
    # Reading shapefiles using GeoPandas. Can take several seconds
    print("@INFOs-----\nReading shapefiles using GeoPandas.")
    # Tested other file format than shp to read vector data
    # it takes 24.38" to read the file in shp format.
    #          23.96"                     gpkg format.
    #          23.96"                     parquet format.
    # -> no significant difference in time to read the file.
    # swtich read_file engine to pyogrio to read the file in 0.9"....

    tiles_df = gpd.read_file(tiles_fn, engine="pyogrio")
    aoi_df = gpd.read_file(args.aoi, engine="pyogrio")

    for i in aoi_df.index:
        print(
            f"@INFOs-----\nIterating through {i+1}/{len(aoi_df)} features within shapefile"
        )
        aoi_row = aoi_df.loc[[i]]
        # Spatial request to select the intersection between two shapefiles
        selection = tiles_df[tiles_df.intersects(aoi_row.geometry.values[0])]
        print(
            f"@INFOs-----\n{len(selection)} tiles intersects '{aoi_row.aoi_name.values[0]}'"
        )

        aoi_path = workdir.joinpath(aoi_row.aoi_name.values[0])
        if not aoi_path.exists():
            # Creating directory if not exist
            aoi_path.mkdir()

        #  Multi-threaded downloading process:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = executor.map(
                lambda i: download_file(
                    selection["url_telech"].values[i], str(extraction_path)
                ),
                range(len(selection)),
            )
            for future in futures:
                try:
                    future.result()  # You can use this to check for download completion or errors
                except FileNotFoundError as exc:
                    print("%r generated a FileNotFoundError:" % exc)
                except Exception as exc:
                    print("%r generated an exception: " % exc)

        # Multi-processing PDAL processing
        cpu_count = psutil.cpu_count()
        max_workers = int(cpu_count * CPU_WORKLOAD)

        list_tiff_files = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            jobs = []
            for j in tqdm(range(len(selection))):
                laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
                laz_fn = (laz_path.name).split(".")[0]
                dem_out_path = aoi_path.joinpath(f"{laz_fn}.tif")
                list_tiff_files.append(str(dem_out_path))
                print(f"Converting '{laz_fn}' file into '{dem_out_path.name}'\n-----")
                job = executor.submit(
                    process_tile, laz_path, dem_out_path, args.dem_resolution
                )
                jobs.append(job)


            for job in concurrent.futures.as_completed(jobs):
                try:
                    job.result()  # Check for processing completion or errors
                except Exception as exc:
                    print(f"Processing generated an exception: {exc}")

        # END for
        # Merging all raster files info single one.
        # TODO: solve problem with border when mosaic tiles.
        os.chdir(aoi_path)
        if args.file_data_type == "gtiff":
            cmd = (
                f"gdal_merge.py -of GTiff -ot Float32 "
                f"-ps {args.dem_resolution} {args.dem_resolution} -n -9999 -a_nodata -9999 "
                f"-o {aoi_df.loc[[i]].aoi_name.values[0]}_{args.dem_resolution}_merged.tif {' '.join(list_tiff_files)}"
            )

        if args.file_data_type == "vrt":
            cmd = (
                f"gdalbuildvrt -tr {args.dem_resolution} {args.dem_resolution} "
                f"-r bilinear {aoi_df.loc[[i]].aoi_name.values[0]}_{args.dem_resolution}_merged.vrt {' '.join(list_tiff_files)}"
            )

        subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        # os.system(cmd)


if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(
        f"Elapsed time in %H:%M:%S: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}"
    )
