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
import sys
import time
import argparse
import subprocess
import shutil
import logging

from pathlib import Path
import wget
import requests
import json
import pdal
from tqdm import tqdm

import geopandas as gpd
import psutil
import concurrent.futures
from shutil import which as find_executable


def getparser():
    args_desc = """
    --------------------------------------------------
    lidarhd downloader --> 'lidarhd_downloader.py' module

    Script to download lidar-hd classified tiles from IGN.
    --------------------------------------------------"""

    # Create an argument parser to get the CLI user arguments
    parser = argparse.ArgumentParser(
        prog="lidarhd_downloader.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=args_desc,
    )
    # MANDATORY - Inputs
    parser.add_argument(
        "code_municipality",
        type=str,
        help="Insee code of the municipality to download LiDAR data.",
    )
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
        "-compute_elev",
        "--compute_elevation",
        default="mean",
        choices=["mean", "min", "max", "median"],
        help="Compute elevation statistics. Default value : mean",
    )
    parser.add_argument(
        "-dtype",
        "--file_data_type",
        type=str,
        default="gtiff",
        choices=["gtif", "vrt"],
        help="Outout data format between GeoTIff and Virtual Dataset (VRT). Default value : gtif",
    )
    parser.add_argument(
        "-force_database",
        "--force_redownload_database",
        action="store_true",
        help="Force re-download of IGN database. Default value: False",
    )
    parser.add_argument(
        "-rm_tiles",
        "--remove_tiles",
        action="store_true",
        help="Remove individual downloaded tiles after processing. Default value: True",
    )
    parser.add_argument(
        "-pdensity",
        "--point_density_map",
        action="store_true",
        help="Generates point density map for given resolution. Requires pdal_wrench installed. Default value: False",
    )

    parser.add_argument(
        "-tin",
        "--triangulation_interpolation",
        action="store_true",
        help="Generates TIN interpolation map for given resolution. Requires pdal_wrench installed. Default value: False",
    )

    parser.add_argument(
        "-cpu_w",
        "--cpu_workload",
        type=float,
        default=0.6,
        help="Multi-threaded process ratio. Default value: 0.6. Max value = 1.0.",
    )
    return parser


# END def


def print_infoBM(text: str, bold: bool = False) -> None:
    """Prints a formatted message in the console.

    Args:
        text (str): string to be printed.
        bold (bool, optional): print the text using bold style. Defaults to False.
    """
    GMC_TEXT = "[ GMC-info ] :"
    if bold:
        print(f"\033[1m{GMC_TEXT} {text}\033[1m")
    else:
        print(f"{GMC_TEXT} {text}")


# END def


def get_municipality_info(code_commune: int):
    """Get the information of a municipality from the IGN API."""
    bbox_cmd = f'curl "https://geo.api.gouv.fr/communes?code={code_commune}&format=geojson&geometry=bbox&fields=code,nom"'
    process = subprocess.Popen(
        bbox_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True
    )
    print(bbox_cmd)
    std_out, std_err = process.communicate()
    if process.returncode != 0:
        print(f"An error occurred: {std_err}")
        return gpd.GeoDataFrame()
    else:
        data = json.loads(std_out)
        gdf = gpd.GeoDataFrame.from_features(data, crs="EPSG:4326")
        gdf.to_crs(epsg=2154, inplace=True)  # Convert the CRS to RGF93 / Lambert-93
        return gdf


def get_municipality_outline(code_commune: int):
    """Get the outline of a municipality from the IGN API."""
    bbox_cmd = f'curl "https://geo.api.gouv.fr/communes?code={code_commune}&format=geojson&geometry=contour&fields=code,nom"'
    print(bbox_cmd)
    process = subprocess.Popen(
        bbox_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True
    )
    std_out, std_err = process.communicate()
    if process.returncode != 0:
        print(f"An error occurred: {std_err}")
        return gpd.GeoDataFrame()
    else:
        data = json.loads(std_out)
        gdf = gpd.GeoDataFrame.from_features(data, crs="EPSG:4326")
        gdf.to_crs(epsg=2154, inplace=True)  # Convert the CRS to RGF93 / Lambert-93
        return gdf


def pdalwrench_bin(bin_name: str) -> bool:
    """Search for PDAL wrench binaries.

    Args:
        bin_name (str): name of binaries to search for.

    Returns:
        bool: boolean value if binaries are found.
    """
    bin_path = find_executable(bin_name)
    if bin_path is None:
        msg = (
            f"Unable to find executable {bin_name}\n"
            f"Install PDAL_WRENCH and ensure it is in your PATH env variable\n"
            "https://github.com/PDAL/wrench"
        )
        sys.exit(msg)

    call = [bin_path]
    print(" ".join(call))

    try:
        result = subprocess.run(
            call, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode == 0:
            return True
        else:
            # print(f"Error: {result.stderr}")
            return False
    except Exception as e:
        print(f"Exception occurred while checking {bin_name}: {e}")
        return False


def pdal_json_pipeline(
    input_laz_fn: str | Path,
    out_tif_fn: str | Path,
    compute_elev: str = "mean",
    tr: float = 1.0,
) -> dict:
    """Generates a PDAL JSON pipeline for processing LiDAR data.

    Args:
        input_laz_fn (str | Path): input laz filename path.
        out_tif_fn (str | Path): Output tif filename path.
        compute_elev (str, optional): parameter for computer elevation. Defaults to "mean". choices are "mean", "min", "max", "median".
        tr (float, optional): resolution of output tif file. Defaults to 1.0.

    Returns:
        dict: PDAL JSON pipeline dictionary.
    """
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


def process_tile(
    laz_path: str | Path,
    output_path: str | Path,
    compute_elev: str = "mean",
    resolution: int | float = 1.0,
) -> None:
    """Process LiDAR tile using PDAL.

    Args:
        laz_path (str | Path): input laz filename path.
        output_path (str | Path): output tif filename path.
        compute_elev (str, optional): parameter for computer elevation. Defaults to "mean". choices are "mean", "min", "max", "median".
        resolution (int | float, optional): resolution of output tif file. Defaults to 1.0.
    """
    # This wraps the PDAL processing
    json_pipeline = pdal_json_pipeline(
        str(laz_path), str(output_path), compute_elev=compute_elev, tr=resolution
    )
    pdal_json_str = json.dumps(json_pipeline)
    pipeline = pdal.Pipeline(pdal_json_str)
    pipeline.execute()
    # END def


def pdal_wrench_density(
    input_laz_fn: str | Path, out_tif: str | Path, tr: int | float = 1.0
) -> list[str]:
    """Compute point density using PDAL Wrench.

    Args:
        input_laz_fn (str | Path): input laz filename path.
        out_tif (str | Path): output tif filename path.
        tr (int | float, optional): resolution of output tif file. Defaults to 1.0.

    Returns:
        list[str]: pdal_wrench command list.
    """
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
    subprocess.run(
        " ".join(cmd_pdal_wrench),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    return cmd_pdal_wrench


def pdal_wrench_tin(
    input_laz_fn: str | Path, out_tif: str | Path, tr: int | float = 1.0
) -> list[str]:
    """Exports filtered point cloud data (class 2 = ground points) to a 2D raster grid using
    a triangulation of points and then interpolating cell values
    from triangles. in short, DEM generation using TIN interpolation.

    Args:
        input_laz_fn (str | Path): input laz filename path.
        out_tif (str | Path): output tif filename path.
        tr (int | float, optional): resolution of output tif file. Defaults to 1.0.

    Returns:
        list[str]: pdal_wrench command list.
    """
    cmd_pdal_wrench = []
    cmd_pdal_wrench.extend(
        [
            "pdal_wrench",
            "to_raster_tin",
            f"--input={input_laz_fn}",
            f"--resolution={tr}",
            f"--filter=Classification==2",
            f"--output={out_tif}",
        ]
    )
    subprocess.run(
        " ".join(cmd_pdal_wrench),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    return cmd_pdal_wrench


def download_file(url: str, output_path: str | Path) -> None:
    """Download file from a given URL.

    Args:
        url (str): URL of the file to download.
        output_path (str | Path): Output path to store the downloaded file.
    """
    print("Start download of ", url)
    response = requests.get(url, timeout=10)
    if "content-disposition" in response.headers:
        content_disposition = response.headers["content-disposition"]
        filename = content_disposition.split("filename=")[1]
        print(filename)
    else:
        filename = url.split("/")[-1]
    if not os.path.exists(os.path.join(output_path, filename)):
        print("Downloading ", filename)
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
        out_laz_dir (str | Pa # Check if directory exist
    # END ifth): Out directory path.
    """
    # Check if file already exist
    if out_dir.joinpath(selected_tiles["nom_pkk"].values[0]).exists():
        print_infoBM(
            f"File {out_dir.joinpath(selected_tiles['nom_pkk'].values[0])} already exist"
        )
        pass
    else:
        print_infoBM(f"Downloading {selected_tiles['nom_pkk'].values[0]}\n-----")
        wget.download(url=selected_tiles["url_telech"].values[0], out=str(out_dir))
    # END if


# END def


def download_LiDAR_tiles_database(out_dir: str | Path) -> None:
    """Download IGN LiDAR database.

    This function will download IGN LiDAR database if not exist.

    Args:
        out_dir (str | Path): Out directory path.
    """
    # First try: Downloading IGN database from official website
    print_infoBM("Checking IGN database from official website . . .")
    try:
        primary_url = (
            "https://diffusion-lidarhd-classe.ign.fr/download/lidar/shp/classe"
        )
        wget.download(url=primary_url, out=str(out_dir))
    except Exception as e:
        print_infoBM(f"Primary download failed with error: {e}")
        print_infoBM("Switching to alternative download link on Zenodo . . .")
        try:
            alternative_url = "https://zenodo.org/records/13793544/files/grille.zip"
            wget.download(url=alternative_url, out=str(out_dir))
        except Exception as e:
            print_infoBM(f"Alternative download failed with error: {e}")
            print_infoBM(
                "Both download attempts failed. Please check your internet connection "
                "or the availability of the URLs."
            )
            exit(1)
    # END try


# END def


def main(args: argparse.Namespace = None):
    if args is None:
        # Get arguments
        parser = getparser()
        args = parser.parse_args()

    project_path = Path(os.environ.get("FIREACCESS_PATH"))
    data_path = project_path / "data"
    data_path.mkdir(parents=True, exist_ok=True)

    raw_lidar_data_path = data_path / "external" / "lidarhd_ign" / "raw_laz_data"
    raw_lidar_data_path.mkdir(parents=True, exist_ok=True)
    lidar_products = data_path / "external" / "lidarhd_ign" / "lidar_products"
    lidar_products.mkdir(parents=True, exist_ok=True)

    # percentage of the CPU workload to be used for processing.
    # Warning: keep some CPUs for the OS and other processes(at least 4 CPUs).
    CPU_WORKLOAD = args.cpu_workload

    if args.out_data_path == None:
        # Default workdir to script's parent directory if out_data_path is not specified
        workdir = data_path / "external" / "lidarhd_ign"
        print_infoBM(
            "As not out_path have been specified, data will be stored by default in 'data/external/lidarhd_ign' directory."
        )
    else:
        workdir = Path(args.out_data_path).resolve()
        if not workdir.exists():
            workdir.mkdir(parents=True)
        # END if
        print_infoBM(f"Data will be stored in {workdir}/raw_laz_data")
    # END if

    # First, we need to download the LiDAR-HD IGN database cmd_pdal_wrench
    print_infoBM("Stage 1 -> Downloading IGN database . . .")

    lidar_database_path = workdir
    print(lidar_database_path)
    if not lidar_database_path.joinpath("TA_programme_LiDAR-HD").exists():
        lidar_database_path.joinpath("TA_programme_LiDAR-HD").mkdir()
    # END if
    tiles_fn = lidar_database_path.joinpath(
        "TA_programme_LiDAR-HD", "TA_diff_pkk_lidarhd_classe.shp"
    )
    if not tiles_fn.exists():
        download_LiDAR_tiles_database(
            lidar_database_path.joinpath("TA_programme_LiDAR-HD")
        )
        # Unzipping downloaded file
        print_infoBM("Unzipping downloaded file . . .")
        zip_database_fn = lidar_database_path.joinpath(
            "TA_programme_LiDAR-HD", "grille.zip"
        )
        result = subprocess.run(
            f"unzip {str(zip_database_fn)} -d {str(zip_database_fn.parent)}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        if result.returncode == 0:
            print_infoBM("Unzipping completed successfully.")
        else:
            print_infoBM(f"Unzipping failed with error: {result.stderr}")
    # END if
    if tiles_fn.exists() & args.force_redownload_database:
        print_infoBM("Forcing to re-download IGN database.")
        ign_ressources_path = workdir.joinpath("TA_programme_LiDAR-HD")
        files_to_delete = ign_ressources_path.glob("*")
        for item in files_to_delete:
            if item.is_dir():
                shutil.rmtree(item)  # Removes directories and their contents
            else:
                item.unlink()  # Removes files
            # END if
        # END for
        download_LiDAR_tiles_database(
            lidar_database_path.joinpath("TA_programme_LiDAR-HD")
        )
    # END if

    print_infoBM(f"Stage 2 -> Working on '{workdir}' directory.")

    # Reading shapefiles using GeoPandas. Can take several seconds
    print_infoBM(f"Reading LiDAR-HD database on {lidar_database_path} . . .")
    tiles_df = gpd.read_file(tiles_fn, engine="pyogrio")

    print_infoBM("Reading AOI info from code commune. . .")

    # # Check if code_com file exists
    code_com = args.code_municipality

    # aoi_df = gpd.read_file(args.aoi_file, engine="pyogrio")
    municipality_gdf = get_municipality_info(code_com)
    roi_name = municipality_gdf.nom.values[0]
    print_infoBM(f"Processing {roi_name} . . .")

    interim_data_path = data_path / "interim" / roi_name
    interim_data_path.mkdir(parents=True, exist_ok=True)
    processed_data_path = data_path / "processed" / roi_name
    processed_data_path.mkdir(parents=True, exist_ok=True)
    extraction_path = workdir.joinpath("raw_laz_data")
    extraction_path.mkdir(parents=True, exist_ok=True)

    municipality_outline_gdf = get_municipality_outline(code_com)
    municipality_outline_gdf.to_file(processed_data_path / f"{roi_name}_contour.gpkg")

    selection = gpd.sjoin(tiles_df, municipality_gdf, how="inner")
    # lst_tiles_id = gdf_tiles_id.nom_left.to_list()

    # for i in aoi_df.index:
    # print_infoBM(f"Iterating through {i+1}/{len(aoi_df)} features within AOI")
    # aoi_row = aoi_df.loc[[i]]
    # # Spatial request to select the intersection between two shapefiles
    # selection = tiles_df[tiles_df.intersects(aoi_row.geometry.values[0])]
    # print_infoBM(
    #     f"{len(selection)} tiles intersects '{aoi_df.loc[[i]].aoi_name.values[0]}' AOI."
    # )
    # aoi_path = interim_data_path

    # END if
    #  Multi-threaded downloading process:
    with concurrent.futures.ThreadPoolExecutor() as executor:
        try:
            for result in executor.map(
                lambda i: download_file(
                    selection["url_telech"].values[i], str(extraction_path)
                ),
                range(len(selection)),
            ):
                # Process result here if needed
                pass
        except FileNotFoundError as exc:
            print_infoBM("%r generated a FileNotFoundError:" % exc)
        except Exception as exc:
            print_infoBM("%r generated an exception: " % exc)

    # Multi-processing PDAL processing
    cpu_count = len(
        psutil.Process().cpu_affinity()
    )  # https://stackoverflow.com/questions/57260410/python-psutil-cpu-count-returns-wrong-number-of-cpus
    max_workers = int(cpu_count * CPU_WORKLOAD)

    print_infoBM("Stage 3 -> Converting *.laz tiles to *.tif DEMs.")
    # Processing LiDAR tiles using multi-threading strategy
    list_tiff_files_merge = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        laz_tif_elevation_jobs = []
        for j in tqdm(range(len(selection))):
            laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
            laz_fn = (laz_path.name).split(".")[0]
            dem_out_path = lidar_products.joinpath(f"{laz_fn}.tif")
            list_tiff_files_merge.append(str(dem_out_path))
            if dem_out_path.exists():
                print_infoBM(
                    f"Skipping conversion for '{laz_path.name}' as output '{dem_out_path.name}' already exists."
                )
                continue

            print_infoBM(
                f"Converting '{laz_path.name}' file into '{dem_out_path.name}'"
            )
            job = executor.submit(
                process_tile,
                laz_path,
                dem_out_path,
                args.compute_elevation,
                args.dem_resolution,
            )
            laz_tif_elevation_jobs.append(job)
        # END for

        for job_a in concurrent.futures.as_completed(laz_tif_elevation_jobs):
            try:
                job_a.result()  # Check for processing completion or errors
            except Exception as exc:
                print_infoBM(f"Processing generated an exception: {exc}")
            # END try
        # END for
    # END with
    # Processing LiDAR point density tiles using multi-threading strategy
    if args.point_density_map:
        print_infoBM("Stage 4 -> Computing LiDAR point density . . .")
        # if pdalwrench_bin("pdal_wrench") is True:
        list_tiff_density_files_merge = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            laz_tif_density_jobs = []
            for j in tqdm(range(len(selection))):
                laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
                laz_fn = (laz_path.name).split(".")[0]
                dem_density_out_path = lidar_products.joinpath(
                    f"{laz_fn}_PointDensity.tif"
                )
                list_tiff_density_files_merge.append(str(dem_density_out_path))
                print_infoBM(
                    f"Computing point density for '{laz_path.name}' file into '{dem_density_out_path.name}'"
                )
                job_density = executor.submit(
                    pdal_wrench_density,
                    str(laz_path),
                    str(dem_density_out_path),
                    args.dem_resolution,
                )
                laz_tif_density_jobs.append(job_density)
            # END for
            for job_b in concurrent.futures.as_completed(laz_tif_density_jobs):
                try:
                    job_b.result()  # Check for processing completion or errors
                except Exception as exc:
                    print_infoBM(f"Processing generated an exception: {exc}")
    # END if
    # Processing LiDAR TIN interpolation tiles using multi-threading strategy
    if args.triangulation_interpolation:
        print_infoBM("Stage 4bis -> Computing LiDAR TIN interpolation . . .")
        list_tiff_tin_files_merge = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            laz_tif_tin_jobs = []
            for j in tqdm(range(len(selection))):
                laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
                laz_fn = (laz_path.name).split(".")[0]
                dem_tin_out_path = lidar_products.joinpath(f"{laz_fn}_TIN.tif")
                list_tiff_tin_files_merge.append(str(dem_tin_out_path))
                # Check if the output file already exists
                if dem_tin_out_path.exists():
                    print_infoBM(
                        f"Skipping TIN interpolation for '{laz_path.name}' as output '{dem_tin_out_path.name}' already exists."
                    )
                    continue

                print_infoBM(
                    f"Computing TIN interpolation for '{laz_path.name}' file into '{dem_tin_out_path.name}'"
                )
                job_tin = executor.submit(
                    pdal_wrench_tin,
                    str(laz_path),
                    str(dem_tin_out_path),
                    args.dem_resolution,
                )
                laz_tif_tin_jobs.append(job_tin)
            for job_c in concurrent.futures.as_completed(laz_tif_tin_jobs):
                try:
                    job_c.result()  # Check for processing completion or errors
                except Exception as exc:
                    print_infoBM(f"Processing generated an exception: {exc}")

    # os.chdir(aoi_path)
    # Merge all tiles by a given resolution
    cmd = []
    merge_out_path = interim_data_path.joinpath(
        f"{roi_name}_Res-{args.dem_resolution}_CompElev-{args.compute_elevation}_merged.tif"
    )
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
                f"{' '.join(list_tiff_files_merge)}",
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
                f"{roi_name}_{args.dem_resolution}_merged.vrt",
                f"{' '.join(list_tiff_files_merge)}",
            ]
        )
    # END if
    print_infoBM(f"Stage 5 -> Merging DEM tiles into '{merge_out_path.name}'.")
    subprocess.run(
        " ".join(cmd),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    # Merge all density tiles by a given resolution
    if args.point_density_map:
        cmd_merge_pointdensity = []
        merge_pointdensity_out_path = interim_data_path.joinpath(
            f"{roi_name}_Res-{args.dem_resolution}_CompElev-{args.compute_elevation}_merged_PointDensity.tif"
        )
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
                f"{' '.join(list_tiff_density_files_merge)}",
            ]
        )
        print_infoBM(
            f"Stage 6 -> Merging LiDAR point density tiles into '{merge_pointdensity_out_path.name}'."
        )
        subprocess.run(
            " ".join(cmd_merge_pointdensity),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

    # Merge all tin tiles by a given resolution
    if args.triangulation_interpolation:
        cmd_merge_tin = []
        merge_tin_out_path = interim_data_path.joinpath(
            f"{roi_name}_Res-{args.dem_resolution}_CompElev-{args.compute_elevation}_merged_TIN.tif"
        )
        cmd_merge_tin.extend(
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
                f"-o  {merge_tin_out_path}",
                f"{' '.join(list_tiff_tin_files_merge)}",
            ]
        )
        print_infoBM(
            f"Stage7 -> Merging LiDAR TIN tiles into '{merge_tin_out_path.name}'."
        )
        subprocess.run(
            " ".join(cmd_merge_tin),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    # TODO: remove tin tiles after processing
    # Remove individual tiles
    if args.remove_tiles:
        print_infoBM("Stage 6 -> Removing individual DEM tiles.")
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
                # END try
            # END for
        # END if
    # END if


# END for


if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(
        f"Elapsed time in %H:%M:%S: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}"
    )
# END if
