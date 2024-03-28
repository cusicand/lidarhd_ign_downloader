#!/usr/bin/env python
# coding=utf-8
#%%
import os
import time
import argparse
import subprocess

from pathlib import Path
import wget
import json
import pdal
from tqdm import tqdm

import geopandas as gpd
#%%
def pdal_json_pipeline(input_laz_fn, out_tif_fn, out_type="mean", tr=1):
    pdal_json_pipeline = {
        "pipeline":
        [
            # To do: check if this is a text or posix object
            input_laz_fn,
        {
            "type":"filters.range",
            "limits":"Classification[1:2]"
        },
        {
            "filename": out_tif_fn,
            "gdaldriver":"GTiff",
            "output_type":out_type,
            "resolution":tr,
            "type": "writers.gdal"
        }
        ]
    }
    return pdal_json_pipeline


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


def download_data(selected_tiles: gpd.GeoDataFrame,
                  out_dir: Path) -> None:
    """Download desired tiles.
    This function check first if data exist into indicated `out_laz_dir`.
    If yes, the function pass the downloading process. If not, will download.
    Args:
        selected_tiles (gpd.GeoDataFrame): Desired ign tiles.
        out_laz_dir (Path): Out directory path.S
    """
    # Check if file already exist
    if out_dir.joinpath(selected_tiles['nom_pkk'].values[0]).exists():
        print(f"File {out_dir.joinpath(selected_tiles['nom_pkk'].values[0])} already exist\n-----")
        pass
    else:
        print(f"Downloading {selected_tiles['nom_pkk'].values[0]}\n-----")
        wget.download(url=selected_tiles['url_telech'].values[0], out=str(out_dir))


def getparser():
    args_desc = """
    --------------------------------------------------
    lidarhd downloader --> 'lidarhd_downloader.py' module

    Script to download lidar-hd classified tiles from IGN.
    --------------------------------------------------"""

    # Create an argument parser to get the CLI user arguments
    parser = argparse.ArgumentParser(
        prog='lidarhd_downloader.py.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=args_desc)
    # MANDATORY - Inputs
    parser.add_argument("aoi", type=str,
        help="path to image aoi.shp file")
    # OPTIONNALS
    parser.add_argument("-out_data", "--out_data_path", type=str, default=None,
        help="Out data-path directory. If not specified, data will be stored in lidar_downloader.py base-path by default.")
    parser.add_argument("-tr", "--dem_resolution", type=float, default=1.0,
        help="resolution of output DEM. Default value : 1.0")
    parser.add_argument("-dtype", "--file_data_type", type=str, default='gtiff', choices=['gtiff', 'vrt'],
        help="Outout data format between GeoTIff and Virtual Dataset (VRT). Default value : GTiff")    
    return parser


def main():
    # Get arguments
    parser = getparser()
    args = parser.parse_args()
    if args.out_data_path == None:
        # workdir = Path.cwd()
        workdir = Path(__file__).parent
        print(f"As not out_path have been specified, data will be stored in $Home/lidarhd_ign_downloader/raw_laz_data\n-----")
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
        wget.download(url="https://diffusion-lidarhd-classe.ign.fr/download/lidar/shp/classe", out=str(workdir.joinpath("ign_resources")))
        subprocess.run(f"unzip {str(workdir.joinpath('ign_resources', 'grille.zip'))} -d {str(workdir.joinpath('ign_resources'))}",
                       shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # os.system(f"unzip {str(workdir.joinpath('ign_resources', 'grille.zip'))} -d {str(workdir.joinpath('ign_resources'))}")

    # Reading shapefiles using GeoPandas. Can take several seconds
    tiles_df = gpd.read_file(tiles_fn)
    aoi_df = gpd.read_file(args.aoi)

    for i in aoi_df.index:
        print(f"Iterating through {i+1}/{len(aoi_df)} features within shapefile\n-----")
        # Spatial request to select the intersection between two shapefiles
        selection = get_lidar_tiles(tiles_df, aoi_df.loc[[i]])
        print(f"{len(selection)} tiles intersects '{aoi_df.loc[[i]].aoi_name.values[0]}'\n-----")
        aoi_path = workdir.joinpath(aoi_df.loc[[i]].aoi_name.values[0])
        # Check if directory exist
        if not aoi_path.exists():
            aoi_path.mkdir()
    
        print(f"{len(selection)} tiles found . . .\n-----")
        for g in range(len(selection)):
            download_data(selection.iloc[[g]], extraction_path)
        # END for
        list_tiff_files = []
        for j in tqdm(range(len(selection))):
            laz_path = extraction_path.joinpath(selection['nom_pkk'].values[j])
            laz_fn = (laz_path.name).split('.')[0]
            dem_out_path = aoi_path.joinpath(f"{laz_fn}.tif")
            print(f"Converting '{laz_fn}' file into '{dem_out_path.name}'\n-----")
            json_pipeline = pdal_json_pipeline(str(laz_path), str(dem_out_path), out_type="mean", tr=args.dem_resolution)
            pdal_json_str = json.dumps(json_pipeline)
            pipeline = pdal.Pipeline(pdal_json_str).execute()
            list_tiff_files.append(str(dem_out_path))

        # Merging all raster files info single one.
        # TODO: solve problem with border when mosaic tiles.
        os.chdir(aoi_path)
        # Merge all tiles by a given resolution

        if args.file_data_type == 'gtiff':
            cmd = f"gdal_merge.py -of GTiff -ot Float32 \
                -ps {args.dem_resolution} {args.dem_resolution} -n -9999 -a_nodata -9999 \
                -o {aoi_df.loc[[i]].aoi_name.values[0]}_{args.dem_resolution}_merged.tif {' '.join(list_tiff_files)}"

        if args.file_data_type == 'vrt':
            cmd = f"gdalbuildvrt -tr {args.dem_resolution} {args.dem_resolution} \
                -r bilinear {aoi_df.loc[[i]].aoi_name.values[0]}_{args.dem_resolution}_merged.vrt {' '.join(list_tiff_files)}"

        subprocess.run(cmd,shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # os.system(cmd)
    #END for

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(f"Elapsed time in %H:%M:%S: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}")
