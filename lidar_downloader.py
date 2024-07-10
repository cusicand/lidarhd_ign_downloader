#!/usr/bin/env python
# coding=utf-8

# ---------------------------------------------------------------------- #
# GeoMultiCorr (GMC) project
# lidar_downloader_jupyter.ipynb
# creation date: 2024-07-09.
#
# Author(s) metadata
# -> author: Diego CUSICANQUI 1 & Jean Baptiste BARRE 2
# -> affiliation 1: CNES | ISTerre | Univ. Grenoble Alpes
# -> affiliation 2: IGE | Univ. Grenoble Alpes
# -> email(s): diego.cusicanqui@univ-grenoble-alpes.fr | diego.cusicanqui.vg@gmail.com
# ->           jb.barre@gmail.com
#
#  Copyright (C) Diego Cusicanqui, 2024 | All rights reserved.
# ---------------------------------------------------------------------- #
# %%
import os
import time
import argparse
import subprocess
import shutil

from pathlib import Path
import wget
import requests
import json
import pdal
from tqdm import tqdm

import geopandas as gpd
import psutil
import concurrent.futures
from distutils.spawn import find_executable

# percentage of the CPU workload to be used for processing. 
# Warning: keep some CPUs for the OS and other processes(at least 4 CPUs).
CPU_WORKLOAD = 0.6
# %%
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
        description=args_desc)
    # MANDATORY - Inputs
    parser.add_argument("aoi_file", type=str, help="path to AOI file. File can be either in *.shp or *.gpkg format")
    # OPTIONNALS
    parser.add_argument("-out_data", "--out_data_path", type=str, default=None,
                        help="Out data-path directory. If not specified, data will be stored in lidar_downloader.py base-path by default.")
    parser.add_argument("-tr", "--dem_resolution", type=float, default=1.0,
                        help="resolution of output DEM. Default value : 1.0")
    parser.add_argument("-compute_elev", "--compute_elevation", default="mean",
                        choices=["mean", "min", "max", "median"],
                        help="Compute elevation statistics. Default value : mean")
    parser.add_argument("-dtype", "--file_data_type", type=str, default="gtiff",
                        choices=["gtif", "vrt"],
                        help="Outout data format between GeoTIff and Virtual Dataset (VRT). Default value : gTif")
    parser.add_argument("-force_database", "--force_redownload_database", action="store_true",
                        help="Force re-download of IGN database. Default value: False")
    parser.add_argument("-rm_tiles", "--remove_tiles", action="store_true",
                        help="Remove individual downloaded tiles after processing. Default value: True")
    parser.add_argument("-pdensity", "--point_density_map", action="store_true",
                        help="Generates point density map for given resolution. Requires pdal_wrench installed. Default value: True")
    return parser
# END def

def print_infoBM(text: str,
                bold: bool = False) -> None:
    GMC_TEXT = "[ GMC-info ] :"
    if bold:
        print(f"\033[1m{GMC_TEXT} {text}\033[1m")
    else:
        print(f"{GMC_TEXT} {text}")
#END def

def pdalwrench_bin(bin, args, **kw) -> bool:
    """Search for PDAL wrench binaries.
    """
    bin_path = find_executable(bin)
    if bin_path is None:
        msg = (f"Unable to find executable {bin_path}\n" 
        f"Install PDAL_WRENCH and ensure it is in your PATH env variable\n" 
        "https://github.com/PDAL/wrench" % bin)
        print(msg)
    try:
        subprocess.run("pdal_wrench", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return True
    except FileNotFoundError:
        return False
# %%
def pdal_json_pipeline(input_laz_fn, out_tif_fn, compute_elev="mean", tr=1.0) -> dict:
    pdal_json_pipeline = {
        "pipeline": [
            # To do: check if this is a text or posix object
            input_laz_fn,
            {"type": "filters.range", "limits": "Classification[0:2]"},
            {
                "filename": out_tif_fn,
                "gdaldriver": "GTiff",
                "output_type": compute_elev,
                "resolution": tr,
                "type": "writers.gdal",
            },
        ]
    }
    return pdal_json_pipeline
# END def

def process_tile(laz_path: str | Path,
                output_path: str | Path,
                compute_elev: str = "mean",
                resolution: int | float = 1.0) -> None:
    # This wraps the PDAL processing
    json_pipeline = pdal_json_pipeline(str(laz_path), str(output_path), compute_elev=compute_elev, tr=resolution)
    pdal_json_str = json.dumps(json_pipeline)
    pipeline = pdal.Pipeline(pdal_json_str)
    pipeline.execute()
    #END def

def pdal_wrench_density(input_laz_fn: str | Path,
                        out_tif: str | Path,
                        tr: int | float = 1.0) -> list[str]:
    cmd_pdal_wrench = []
    cmd_pdal_wrench.extend(
        [
            "pdal_wrench",
            "density",
            f"--input={input_laz_fn}",
            f"--resolution={tr}",
            f"--output={out_tif}",
        ]
    )
    subprocess.run(' '.join(cmd_pdal_wrench), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return cmd_pdal_wrench

def download_file(url: str,
                output_path: str | Path) -> None:
    response = requests.get(url, timeout=10)
    if "content-disposition" in response.headers:
        content_disposition = response.headers["content-disposition"]
        filename = content_disposition.split("filename=")[1]
    else:
        filename = url.split("/")[-1]
    if not os.path.exists(os.path.join(output_path, filename)):
        with open(os.path.join(output_path, filename), mode="wb") as file:
            file.write(response.content)
        print_infoBM(f"{filename} Downloaded.")
    else:
        print_infoBM(f"{filename} already exists.")

def get_lidar_tiles(
    tiles_df: gpd.GeoDataFrame, aoi_feature: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Spatial intersection of two polygons

    Args:
        tiles_df (gpd.GeoDataFrame): LiDAR ign reference tiles.
        aoi_feature (gpd.GeoDataFrame): aoi single feature.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame of intersected features
    """
    # Spatial request to select the intersection between two shapefiles
    intersection = tiles_df[tiles_df.intersects(aoi_feature.geometry.values[0])]
    return intersection
# END def

def download_data(selected_tiles: gpd.GeoDataFrame, out_dir: str | Path) -> None:
    """Download desired tiles.
    This function check first if data exist into indicated `out_laz_dir`.
    If yes, the function pass the downloading process. If not, will download.
    Args:
        selected_tiles (gpd.GeoDataFrame): Desired ign tiles.
        out_laz_dir (Path): Out directory path.S
    """
    # Check if file already exist
    if out_dir.joinpath(selected_tiles["nom_pkk"].values[0]).exists():
        print_infoBM(f"File {out_dir.joinpath(selected_tiles['nom_pkk'].values[0])} already exist")
        pass
    else:
        print_infoBM(f"Downloading {selected_tiles['nom_pkk'].values[0]}\n-----")
        wget.download(url=selected_tiles["url_telech"].values[0], out=str(out_dir))
    # END if
# ENd def

# %%
def main(args: argparse.Namespace = None):
    if args is None:
        # Get arguments
        parser = getparser()
        args = parser.parse_args()
    # END if
    if args.out_data_path == None:
        # Default workdir to script's parent directory if out_data_path is not specified
        workdir = Path(__file__).resolve().parent
        print_infoBM(f"As not out_path have been specified, data will be stored by default in $HOME/lidarhd_ign_downloader/raw_laz_data")
    else:
        workdir = Path(args.out_data_path).resolve()
        if not workdir.exists():
            workdir.mkdir(parents=True)
        # END if
        print_infoBM(f"Data will be stored in {workdir}/raw_laz_data")
    # END if
    print_infoBM(f"Working on: {workdir}")

    extraction_path = workdir.joinpath("raw_laz_data")
    if not extraction_path.exists():
        # Creating directory if not exist
        extraction_path.mkdir(parents=True)
    # END if
    tiles_fn = workdir.joinpath("ign_resources", "TA_diff_pkk_lidarhd_classe.shp")

    if not workdir.joinpath("ign_resources").exists():
        workdir.joinpath("ign_resources").mkdir()
    # END if

    if not tiles_fn.exists():
        print_infoBM("Downloading IGN database . . .")
        wget.download(
            url="https://diffusion-lidarhd-classe.ign.fr/download/lidar/shp/classe",
            out=str(workdir.joinpath("ign_resources")))
        subprocess.run(
            f"unzip {str(workdir.joinpath('ign_resources', 'grille.zip'))} -d {str(workdir.joinpath('ign_resources'))}",
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if tiles_fn.exists() & args.force_redownload_database:
        print_infoBM("Forcing to re-download IGN database.")
        ign_ressources_path = workdir.joinpath("ign_resources")
        files_to_delete = ign_ressources_path.glob("*")
        for item in files_to_delete:
            if item.is_dir():
                shutil.rmtree(item)  # Removes directories and their contents
            else:
                item.unlink()  # Removes files
            #END if
        #END for
        wget.download(
            url="https://diffusion-lidarhd-classe.ign.fr/download/lidar/shp/classe",
            out=str(workdir.joinpath("ign_resources")))
        subprocess.run(
            f"unzip {str(workdir.joinpath('ign_resources', 'grille.zip'))} -d {str(workdir.joinpath('ign_resources'))}",
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # END if

    # Reading shapefiles using GeoPandas. Can take several seconds
    print_infoBM("Reading lidarHD database . . .")
    tiles_df = gpd.read_file(tiles_fn, engine="pyogrio")
    print_infoBM("Reading AOI file . . .")
    aoi_df = gpd.read_file(args.aoi_file, engine="pyogrio")

    for i in aoi_df.index:
        print_infoBM(f"Iterating through {i+1}/{len(aoi_df)} features within AOI")
        aoi_row = aoi_df.loc[[i]]
        # Spatial request to select the intersection between two shapefiles
        selection = tiles_df[tiles_df.intersects(aoi_row.geometry.values[0])]
        print_infoBM(f"{len(selection)} tiles intersects '{aoi_df.loc[[i]].aoi_name.values[0]}' AOI.")
        aoi_path = workdir.joinpath(aoi_df.loc[[i]].aoi_name.values[0])
        # Check if directory exist
        if not aoi_path.exists():
            aoi_path.mkdir()
        # END if

        #  Multi-threaded downloading process:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            try:
                for result in executor.map(
                    lambda i: download_file(selection["url_telech"].values[i], str(extraction_path)),
                    range(len(selection))):
                    # Process result here if needed
                    pass
            except FileNotFoundError as exc:
                print_infoBM("%r generated a FileNotFoundError:" % exc)
            except Exception as exc:
                print_infoBM("%r generated an exception: " % exc)

        # Multi-processing PDAL processing
        cpu_count = psutil.cpu_count()
        max_workers = int(cpu_count * CPU_WORKLOAD)
        
        list_tiff_files_merge = []
        list_tiff_density_files_merge = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            laz_tif_elevation_jobs = []
            laz_tif_density_jobs = []
            for j in tqdm(range(len(selection))):
                laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
                laz_fn = (laz_path.name).split(".")[0]
                dem_out_path = aoi_path.joinpath(f"{laz_fn}.tif")
                dem_density_out_path = aoi_path.joinpath(f"{laz_fn}_PointDensity.tif")
                list_tiff_files_merge.append(str(dem_out_path))
                list_tiff_density_files_merge.append(str(dem_density_out_path))
                print_infoBM(f"Converting '{laz_path.name}' file into '{dem_out_path.name}'")

                job = executor.submit(
                    process_tile, laz_path, dem_out_path, args.compute_elevation, args.dem_resolution)
                laz_tif_elevation_jobs.append(job)
                if args.point_density_map:
                    print_infoBM(f"Computing point density for '{laz_path.name}' file into '{dem_density_out_path.name}'")
                    job_density = executor.submit(
                        pdal_wrench_density, str(laz_path), str(dem_density_out_path), args.dem_resolution)
                    laz_tif_density_jobs.append(job_density)
            # END for
            
            for job_a in concurrent.futures.as_completed(laz_tif_elevation_jobs):
                try:
                    job_a.result()  # Check for processing completion or errors
                except Exception as exc:
                    print_infoBM(f"Processing generated an exception: {exc}")
                #END try
            # END for
            if args.point_density_map:    
                for job_b in concurrent.futures.as_completed(laz_tif_density_jobs):
                    try:
                        job_b.result()  # Check for processing completion or errors
                    except Exception as exc:
                        print_infoBM(f"Processing generated an exception: {exc}")
                    #END try
                # END for
        # TODO: solve problem with border when mosaic tiles.
        os.chdir(aoi_path)
        # Merge all tiles by a given resolution
        cmd = []
        merge_out_path = aoi_path.joinpath(f"{aoi_df.loc[[i]].aoi_name.values[0]}_Res{args.dem_resolution}_CompElev{args.compute_elevation}_merged.tif")
        if args.file_data_type == "gtif":
            cmd.extend(
                [
                    "gdal_merge.py",
                    "-of",
                    "GTiff",
                    "-ot",
                    "Float32",
                    f"-ps {args.dem_resolution} {args.dem_resolution}",
                    "-n",
                    "-9999",
                    "-a_nodata",
                    "-9999",
                    f"-o {merge_out_path}",
                    f"{' '.join(list_tiff_files_merge)}"
                ]
            )
        if args.file_data_type == "vrt":
            cmd.extend(
                [
                    "gdalbuildvrt",
                    "-tr",
                    f"{args.dem_resolution}",
                    f"{args.dem_resolution}",
                    "-r",
                    "bilinear",
                    f"{aoi_df.loc[[i]].aoi_name.values[0]}_{args.dem_resolution}_merged.vrt",
                    f"{' '.join(list_tiff_files_merge)}"
                ]
            )
        #END if
        print_infoBM(f"Merging DEM tiles into '{merge_out_path.name}'.")
        subprocess.run(
            ' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        # Merge all density tiles by a given resolution
        if args.point_density_map:
            cmd_merge_pointdensity = []
            merge_pointdensity_out_path = aoi_path.joinpath(f"{aoi_df.loc[[i]].aoi_name.values[0]}_Res{args.dem_resolution}_CompElev{args.compute_elevation}_merged_PointDensity.tif")
            cmd_merge_pointdensity.extend(
                [
                    "gdal_merge.py",
                    "-of",
                    "GTiff",
                    "-ot",
                    "Float32",
                    f"-ps {args.dem_resolution} {args.dem_resolution}",
                    "-n",
                    "-9999",
                    "-a_nodata",
                    "-9999",
                    f"-o  {merge_pointdensity_out_path}",
                    f"{' '.join(list_tiff_density_files_merge)}"
                ]
            )
            print_infoBM(f"Merging LiDAR point density tiles into '{merge_pointdensity_out_path.name}'.")
            subprocess.run(
                ' '.join(cmd_merge_pointdensity),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
        # Remove individual tiles
        if args.remove_tiles:
            print_infoBM("Removing individual DEM tiles.")
            for dem_tile in list_tiff_files_merge:
                try:
                    os.remove(dem_tile)
                except FileNotFoundError:
                    raise ValueError(f"File {dem_tile} not found. Cannot remove.")
                    pass
            # END for
            if args.point_density_map:
                print_infoBM("Removing individual LiDAR point density tiles.")
                for density_tile in list_tiff_density_files_merge:
                    try:
                        os.remove(density_tile)
                    except FileNotFoundError:
                        raise ValueError(f"File {density_tile} not found. Cannot remove.")
                        pass
            # END for
        # END if
    # END for
# %%
if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(f"Elapsed time in %H:%M:%S: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}")
#END if