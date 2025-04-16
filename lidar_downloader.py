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

import requests
import json
import pdal
from tqdm import tqdm
from pathlib import Path
from datetime import datetime

import pandas as pd
import geopandas as gpd
import psutil
import concurrent.futures
from distutils.spawn import find_executable
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
    # parser.add_argument("aoi_file", type=str, help="path to AOI file. File can be either in *.shp or *.gpkg format")
    parser.add_argument("aoi_file", type=str, nargs='+', help="path(s) to AOI file(s). Files can be either in *.shp or *.gpkg format")
    # OPTIONNALS
    parser.add_argument("-out_data", "--out_data_path", type=str, default=None,
                        help="Out data-path directory. If not specified, data will be stored in lidar_downloader.py base-path by default.")
    parser.add_argument("-tr", "--dem_resolution", type=float, default=1.0,
                        help="resolution of output DEM. Default value : 1.0")
    parser.add_argument("-compute_elev", "--compute_elevation", default="mean",
                        choices=["mean", "min", "max", "median"],
                        help="Compute elevation statistics. Default value : mean")
    parser.add_argument("-dtype", "--file_data_type", type=str, default="gtif",
                        choices=["gtif", "vrt"],
                        help="Outout data format between GeoTIff and Virtual Dataset (VRT). Default value : gtif")
    parser.add_argument("-force_database", "--force_redownload_database", action="store_true",
                        help="Force re-download of IGN database. Default value: False")
    parser.add_argument("-rm_tiles", "--remove_tiles", action="store_true",
                        help="Remove individual downloaded tiles after processing. Default value: True")
    parser.add_argument("-pdensity", "--point_density_map", action="store_true",
                        help="Generates point density map for given resolution. Requires pdal_wrench installed. Default value: False")
    parser.add_argument("-cpu_w", "--cpu_workload", type=float, default=0.6,
                        help="Multi-threaded processing ratio. The maximum value of 0.8 is used to ensure that the programme does not crash. Default value: 0.6. Max value = 0.8.")
    return parser
# END def
#%%
def print_infoBM(text: str,
                bold: bool = False) -> None:
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
#END def
#%%
def pdalwrench_bin(bin_name: str) -> bool:
    """Search for PDAL wrench binaries.

    Args:
        bin_name (str): name of binaries to search for.

    Returns:
        bool: boolean value if binaries are found.
    """
    bin_path = find_executable(bin_name)
    if bin_path is None:
        msg = (f"Unable to find executable {bin_name}\n" 
               f"Install PDAL_WRENCH and ensure it is in your PATH env variable\n" 
               "https://github.com/PDAL/wrench")
        sys.exit(msg)
    
    call = [bin_path]
    print(' '.join(call))
    
    try:
        result = subprocess.run(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            return True
        else:
            # print(f"Error: {result.stderr}")
            return False
    except Exception as e:
        print(f"Exception occurred while checking {bin_name}: {e}")
        return False
# %%
def pdal_json_pipeline(input_laz_fn: str | Path,
                       out_tif_fn: str | Path,
                       compute_elev: str = "mean",
                       tr: float = 1.0) -> dict:
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
            input_laz_fn,
            {"type": "filters.range",
             "limits": "Classification[0:2]"}, #! TODO: check for better classification
            {
                "filename": out_tif_fn,
                "gdaldriver": "GTiff",
                "output_type": compute_elev,
                "resolution": tr,
                "type": "writers.gdal",
                # #! TODO: check for better compression options
                # "creation_options": [
                #     "COMPRESS=LZW",
                #     "PREDICTOR=2",
                #     "ZLEVEL=9",
                #     "NUM_THREADS=ALL_CPUS",
                # ],
                # "nodata": -9999,
                # "data_type": "Float32",
            },
        ]
    }
    return pdal_json_pipeline
# END def
#%%
def process_tile(laz_path: str | Path,
                output_path: str | Path,
                compute_elev: str = "mean",
                resolution: int | float = 1.0) -> None:
    """Process LiDAR tile using PDAL.

    Args:
        laz_path (str | Path): input laz filename path.
        output_path (str | Path): output tif filename path.
        compute_elev (str, optional): parameter for computer elevation. Defaults to "mean". choices are "mean", "min", "max", "median".
        resolution (int | float, optional): resolution of output tif file. Defaults to 1.0.
    """
    # This wraps the PDAL processing
    json_pipeline = pdal_json_pipeline(str(laz_path), str(output_path), compute_elev=compute_elev, tr=resolution)
    pdal_json_str = json.dumps(json_pipeline)
    pipeline = pdal.Pipeline(pdal_json_str)
    pipeline.execute()
    #END def
#%%
def pdal_wrench_density(input_laz_fn: str | Path,
                        out_tif: str | Path,
                        tr: int | float = 1.0) -> list[str]:
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
    subprocess.run(' '.join(cmd_pdal_wrench), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return cmd_pdal_wrench
#%%
def download_file(url: str,
                output_path: str | Path) -> None:
    """Download file from a given URL.

    Args:
        url (str): URL of the file to download.
        output_path (str | Path): Output path to store the downloaded file.
    """
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

def get_lidar_tiles(tiles_df: gpd.GeoDataFrame,
                    aoi_feature: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
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
#%%
def fetch_chunk(url: str,
                ntiles: int,
                start_index: int) -> gpd.GeoDataFrame:
    """Fetch a chunk of data from the given URL.
    This function constructs a URL with the specified parameters and retrieves
    the corresponding data. It handles any exceptions that may occur during
    the retrieval process and returns None if an error occurs.

    Args:
        url (str): tiles url
        ntiles (int): number of tiles to fetch
        start_index (int): index to start fetching from

    Returns:
        gpd.GeoDataFrame: returns a GeoDataFrame containing the fetched data, or None if an error occurs.
    """
    try:
        params = f"&STARTINDEX={start_index}&COUNT={ntiles}&SRSNAME=urn:ogc:def:crs:EPSG::2154"
        full_url = url + params
        gdf = gpd.read_file(full_url)
        return gdf if not gdf.empty else None
    except Exception as e:
        # print(f"WARNING: Error fetching index {start_index}: {e}")
        return None
# END def
# %%
def url2bloc(url_series: pd.Series) -> pd.Series:
    """Convert URLs to block identifiers.
    Args:
        url_series (pd.Series): Series containing URLs.
    Returns:
        pd.Series: Series containing block identifiers.
    """
    # Example dummy implementation
    return url_series.apply(lambda x: x.split("/")[-1].split(".")[0] if isinstance(x, str) else None)
# END def
# %%
def download_lidarhd_database(
    url: str = "https://data.geopf.fr/private/wfs/wfs?apikey=interface_catalogue&SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=IGNF_LIDAR-HD_TA:nuage-dalle",
    out_dir: str | Path = None,
    max_pages: int = 100,
    ntiles: int =5000,
    cpu_workers: int | float = 12) -> str | Path:
    """Downloads the LiDAR-HD database from the specified URL.
    This function uses a ThreadPoolExecutor to fetch data in parallel,
    and it concatenates the fetched data into a single GeoDataFrame.
    It also processes the data by renaming columns and dropping unnecessary ones.
    Finally, it saves the processed data to a GeoPackage file in the specified output directory.
    If the output directory is not specified, it defaults to the current working directory.
    The function returns the path to the saved GeoPackage file.

    Args:
        url (_type_, optional): WFS url. Defaults to "https://data.geopf.fr/private/wfs/wfs?apikey=interface_catalogue&SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=IGNF_LIDAR-HD_TA:nuage-dalle".
        out_dir (str | Path, optional): outout directory. Defaults to None.
        max_pages (int, optional): max number of pages of WFS service. Defaults to 100.
        ntiles (int, optional): number of tiles. Defaults to 5000.
        cpu_workers (int | float, optional): number of workers for parallel processing. Defaults to 1.

    Returns:
        str | Path: returns database filename path.
    """
    out_dir = Path(out_dir)
    start_indexes = [(n - 1) * ntiles for n in range(1, max_pages + 1)]
    data_chunks = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=cpu_workers) as executor:
        futures = {executor.submit(fetch_chunk, url, ntiles, idx): idx for idx in start_indexes}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is not None:
                data_chunks.append(result)
    #END with
    # Concatenate all data chunks into a single DataFrame and process it
    if data_chunks:
        database = pd.concat(data_chunks, ignore_index=True)
        database['bloc'] = url2bloc(database['url'])
        database = database.rename(columns=lambda x: x.replace("url", "url_telech").replace("name", "nom_pkk"))
        database = database.drop(columns=['gml_id'], errors='ignore')
    #END if
    database_filename_path = out_dir / f"LidarHD_tiles_database_{datetime.today().strftime('%Y-%m-%d')}.gpkg"
    print(f"Saving database in: {database_filename_path}")
    database.to_file(out_dir / f"{database_filename_path}", driver='GPKG')
    #END if
    return database_filename_path
#END def
# %%
def main(args: argparse.Namespace = None):
    if args is None:
        # Get arguments
        parser = getparser()
        args = parser.parse_args()
    # END if
    aois_files_list = [Path(aoi_file).resolve() for aoi_file in args.aoi_file]
    # percentage of the CPU workload to be used for processing. 
    # Warning: keep some CPUs for the OS and other processes(at least 4 CPUs).
    if args.cpu_workload > 0.8:
        CPU_WORKLOAD = 0.8
    else:
        CPU_WORKLOAD = args.cpu_workload
    # END if    
    # Multi-processing PDAL processing
    cpu_count = len(psutil.Process().cpu_affinity()) # https://stackoverflow.com/questions/57260410/python-psutil-cpu-count-returns-wrong-number-of-cpus
    max_workers = int(cpu_count * CPU_WORKLOAD)
    
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
    # First, we need to download the LiDAR-HD IGN database
    print_infoBM("Stage 1 -> Looking for LiDAR-HD database . . .")
    lidar_database_path = Path(__file__).resolve().parent
    # Check if ign_resources directory exist
    # print(lidar_database_path)
    if not lidar_database_path.joinpath("ign_resources").exists():
        lidar_database_path.joinpath("ign_resources").mkdir()
    # END if
    tiles_fn = list(lidar_database_path.joinpath("ign_resources").glob("LidarHD_tiles_database_????-??-??.gpkg"))
    if not tiles_fn:
        with tqdm(total=1,
                  desc=f"{print_infoBM('No LiDAR HD database found. Downloading from WFS service . . .')}",
                  unit="task") as pbar:
            tiles_fn = download_lidarhd_database(out_dir=lidar_database_path / "ign_resources")
            pbar.update(1)
    else:
        print_infoBM(f"LiDAR HD database found in {str(tiles_fn)}.")
        tiles_fn = Path(tiles_fn[0])
        pass
        if args.force_redownload_database:
            print_infoBM("Forcing to re-download LiDAR HD database.")
            tiles_fn[0].unlink(missing_ok=True)  # Remove the existing file
            with tqdm(total=1,
                      desc=f"{print_infoBM('Downloading from WFS service . . .')}",
                      unit="task") as pbar:
                tiles_fn = download_lidarhd_database(out_dir=lidar_database_path / "ign_resources")
                pbar.update(1)
        # END if
    # END if
    print_infoBM(f"Stage 2 -> Working on '{workdir}' directory.")

    extraction_path = workdir.joinpath("raw_laz_data")
    if not extraction_path.exists():
        # Creating directory if not exist
        extraction_path.mkdir(parents=True)
    # END if

    # Reading shapefiles using GeoPandas. Can take several seconds
    print_infoBM(f"Reading LiDAR-HD database in {tiles_fn} . . .")
    tiles_df = gpd.read_file(str(tiles_fn), engine="pyogrio")

    for aoi_file in aois_files_list:
        print_infoBM(f"Intersecting tiles for '{aoi_file}' AOI file . . .")
        # print_infoBM(f"Reading {len(aois_files_list)} AOI files . . .")
        aoi_df = gpd.read_file(aoi_file, engine="pyogrio")
        if aoi_df.crs.to_epsg() != 2154:
            # Reprojecting AOI file to match LiDAR-HD database CRS
            print_infoBM(f"Reprojecting AOI file {aoi_file} to match LiDAR-HD database CRS.")
            aoi_df = aoi_df.to_crs(epsg=2154)
        # END if
        selection = tiles_df[tiles_df.intersects(aoi_df.geometry.values[0])]
        print_infoBM(f"{len(selection)} tiles intersects '{aoi_file.name}' AOI.")
        aoi_path = workdir.joinpath(aoi_file.name.split(".")[0])
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
        # END with        
        print_infoBM(f"Stage 3 -> Converting *.laz tiles to *.tif DEMs.")
        # Processing LiDAR tiles using multi-threading strategy
        list_tiff_files_merge = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            laz_tif_elevation_jobs = []
            for j in tqdm(range(len(selection))):
                laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
                laz_fn = (laz_path.name).split(".")[0]
                dem_out_path = aoi_path.joinpath(f"{laz_fn}.tif")
                list_tiff_files_merge.append(str(dem_out_path))
                print_infoBM(f"Converting '{laz_path.name}' file into '{dem_out_path.name}'")
                job = executor.submit(
                    process_tile, laz_path, dem_out_path, args.compute_elevation, args.dem_resolution)
                laz_tif_elevation_jobs.append(job)
            # END for
            
            for job_a in concurrent.futures.as_completed(laz_tif_elevation_jobs):
                try:
                    job_a.result()  # Check for processing completion or errors
                except Exception as exc:
                    print_infoBM(f"Processing generated an exception: {exc}")
                #END try
            # END for
        #END with
        # Processing LiDAR point density tiles using multi-threading strategy
        if args.point_density_map:
            print_infoBM("Stage 4 -> Computing LiDAR point density . . .")
            # if pdalwrench_bin("pdal_wrench") is True:
            list_tiff_density_files_merge = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                laz_tif_density_jobs = []
                for j in tqdm(range(len(selection))):
                    laz_path = extraction_path.joinpath(selection["nom_pkk"].values[j])
                    laz_fn = (laz_path.name).split(".")[0]
                    dem_density_out_path = aoi_path.joinpath(f"{laz_fn}_PointDensity.tif")
                    list_tiff_density_files_merge.append(str(dem_density_out_path))
                    print_infoBM(f"Computing point density for '{laz_path.name}' file into '{dem_density_out_path.name}'")
                    job_density = executor.submit(
                        pdal_wrench_density, str(laz_path), str(dem_density_out_path), args.dem_resolution)
                    laz_tif_density_jobs.append(job_density)
                # END for
                for job_b in concurrent.futures.as_completed(laz_tif_density_jobs):
                    try:
                        job_b.result()  # Check for processing completion or errors
                    except Exception as exc:
                        print_infoBM(f"Processing generated an exception: {exc}")
                    #END try
                # END for
            #END with
        #END if
        os.chdir(aoi_path)
        # Merge all tiles by a given resolution
        cmd = []
        merge_out_path = aoi_path.joinpath(f"{aoi_file.name.split('.')[0]}_LiDARDEM_Res-{args.dem_resolution}_CompElev-{args.compute_elevation}_merged.tif")
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
                    f"{aoi_file.name.split('.')[0]}_{args.dem_resolution}_merged.vrt",
                    f"{' '.join(list_tiff_files_merge)}"
                ]
            )
        #END if
        print_infoBM(f"Stage 5 -> Merging DEM tiles into '{merge_out_path.name}'.")
        subprocess.run(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Merge all density tiles by a given resolution
        if args.point_density_map:
            cmd_merge_pointdensity = []
            merge_pointdensity_out_path = aoi_path.joinpath(f"{aoi_file.name.split('.')[0]}_LiDARDEM_Res-{args.dem_resolution}_CompElev-{args.compute_elevation}_merged_PointDensity.tif")
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
            print_infoBM(f"Stage 6 -> Merging LiDAR point density tiles into '{merge_pointdensity_out_path.name}'.")
            subprocess.run(' '.join(cmd_merge_pointdensity),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
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
                    #END try
                #END for
            # END if
        # END if
    # END for
    print_infoBM("Processing done. Well done bibi!!")
# %%
if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(f"Elapsed time in %H:%M:%S: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}")
#END if