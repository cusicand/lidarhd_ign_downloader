#!/usr/bin/env python
# coding=utf-8
#%%
import os
import sys
import time
import argparse

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
            "limits":"Classification[0:2]"
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
#%%
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
    parser.add_argument("-tr", "--dem_resolution", type=float, default=1.0,
        help="resolution of output DEM. Default value : 1.0")
    return parser
#%%
def main():
    # Get arguments
    parser = getparser()
    args = parser.parse_args()

    # workdir = Path.cwd()
    workdir = Path(__file__).parent
    print(f"@INFOs-----\nWe assume that all data will be stored in $Home/lidarhd_ign_downloader/extractions\nWorking on: {workdir}\n@INFOs-----")

    # We assume that all data will be stored in $Home/lidarhd_ign_downloader/extractions
    # home_path = f"path/to/the/other/location"
    # home_path = workdir.home()
    extraction_path = workdir.joinpath("extractions")
    if not extraction_path.exists():
        # Creating directory if not exist
        extraction_path.mkdir()
    #END if
    tiles_fn = workdir.joinpath("resources", "TA_diff_pkk_lidarhd_classe.shp")

    if not workdir.joinpath("resources").exists():
        workdir.joinpath("resources").mkdir()

    if not tiles_fn.exists():
         wget.download(url="https://diffusion-lidarhd-classe.ign.fr/download/lidar/shp/classe", out=str(workdir.joinpath("resources")))
         os.system(f"unzip {str(workdir.joinpath('resources', 'grille.zip'))} -d {str(workdir.joinpath('resources'))}")
    # Reading shapefiles using GeoPandas. Can take several seconds
    tiles_df = gpd.read_file(tiles_fn)
    aoi_df = gpd.read_file(args.aoi)

    for i in aoi_df.index:
        print(f"@INFOs-----\nIterating through {i+1}/{len(aoi_df)} features within shapefile\n@INFOs-----")
        aoi_row = aoi_df.loc[[i]]
        # Spatial request to select the intersection between two shapefiles
        selection = tiles_df[tiles_df.intersects(aoi_row.geometry.values[0])]
        # selection_list = selection['url_telech'].to_list()
        print(f"@INFOs-----\n{len(selection)} tiles intersects '{aoi_row.aoi_name.values[0]}'\n@INFOs-----")
        # aoi_path = workdir.joinpath(args.outdir)
        aoi_path = workdir.joinpath(aoi_row.aoi_name.values[0])
        if not aoi_path.exists():
            # Creating directory if not exist
            aoi_path.mkdir()
        #END if
        #%%
        for i in range(len(selection)):
            #TO DO: put all the file in the same folder. Ask to more space in the server.
            if extraction_path.joinpath(selection['nom_pkk'].values[i]).exists():
                print(f"File {extraction_path.joinpath(selection['nom_pkk'].values[i])} already exist")
            else:
                print(f"Downloading {i}/{len(selection)} : '{selection['nom_pkk'].values[i]}'")
                wget.download(url=selection['url_telech'].values[i], out=str(extraction_path))
            # END if
        # ENd for
        for j in tqdm(range(len(selection))):
            laz_path = extraction_path.joinpath(selection['nom_pkk'].values[j])
            laz_fn = (laz_path.name).split('.')[0]
            dem_out_path = aoi_path.joinpath(f"{laz_fn}.tif")
            print(f"@INFOs-----\nConverting '{laz_fn}' file into '{dem_out_path.name}'\n@INFOs-----")
            json_pipeline = pdal_json_pipeline(str(laz_path), str(dem_out_path), out_type="mean", tr=args.dem_resolution)
            pdal_json_str = json.dumps(json_pipeline)
            pipeline = pdal.Pipeline(pdal_json_str)
            count = pipeline.execute()
        #END for
        # Merging all raster files info single one.
        os.chdir(aoi_path)
        cmd = f"gdal_merge.py -of GTiff -ot Float32 \
            -ps {args.dem_resolution} {args.dem_resolution} \
            -n -9999 -a_nodata -9999 \
            -o {aoi_row.aoi_name.values[0]}_{args.dem_resolution}_merged.tif *.tif"
        os.system(cmd)
    #END for
# %%
if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(f"Elapsed time in %H:%M:%S: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}")
# %%