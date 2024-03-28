
import os
import sys
import time

from pathlib import Path
from tqdm import tqdm
import geopandas as gpd

shp_file='./ign_resources/TA_diff_pkk_lidarhd_classe.shp'
gkpg_file='./ign_resources/TA_diff_pkk_lidarhd_classe.gpkg'
parquet_file='./ign_resources/TA_diff_pkk_lidarhd_classe.parquet'

def read_vector_file(file_path):
        return gpd.read_file(file_path, engine="pyogrio")
    

def main():
    # time to read shapefile using GeoPandas
    
    start_time = time.time()
    tiles_shp_df = read_vector_file(shp_file)
    elapsed_time = time.time() - start_time 
    print(f"Elapsed time to read shapefile using GeoPandas: {elapsed_time}")
    
    #tiles_shp_df.to_parquet(parquet_file)
    

    

if __name__ == "__main__":
    main()
    
    
