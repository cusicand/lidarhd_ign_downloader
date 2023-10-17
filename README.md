# IGN LiDQR HD downloader
Python script to download and resample LiDAR HD pre-clasified tiles from IGN.

# Required packages
- Python 3.9
- PDAL
- Geopandas
- wget
- pathlib2

If you have a proper python environment already installed, you can install the packages directly by typing:
`conda install -c conda-forge pdal wget pathlib2 geopandas`

# How to use

Type `lidar_downloader.py -h` to access to the help of the tool.

This small tool needs and Area of Interest (AOI) in shapefile format as mandatory argument.
Then, we can switch between several parameters like:
- `-o` or `--outdir` to specify the output directory.
- `-tr` or `--dem_resolution` to select the desired resolution of the output DEM.

# Installation
if you want to create a specific environment, please follow the instructions below

