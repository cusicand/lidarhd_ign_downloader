# IGN LiDAR HD downloader
Python script to download and resample LiDAR HD pre-clasified tiles from IGN.

## Installation

In order to install properly, you have to clone RGDyn repository.
```
git clone https://github.com/cusicand/lidarhd_ign_downloader.git
cd ./lidarhd_ign_downloader
```

## Make python script executables

For both cases, you have to add the following lines to your `.bashrc` file to have complete access to the developed scripts since command line.
```
vi ~/.bashrc or nano ~/.bashrc
export LIDAR_PATH=$HOME/lidarhd_ign_downloader
export PATH=$LIDAR_PATH:$PATH            
export PYTHONPATH=$LIDAR_PATH:$PYTHONPATH

chmod +x lidar_downloader.py
```
Use `source ~/.bashrc` to reload changes.

If you have a proper python environment already installed, you can install the packages directly by typing:
`conda install -c conda-forge pdal wget pathlib2 geopandas`

## How to use

Type `python lidar_downloader.py -h` to access to the help of the tool.

This small tool needs and Area of Interest (AOI) in shapefile (`some_amazing_name.shp`) format as mandatory argument.
Then, we can switch between several parameters like:.
- `-tr` or `--dem_resolution` to select the desired resolution of the output DEM.

below an example of how to extract data by using an AOI file:

`python lidarhd_downloader.py some_amazing_name.shp -tr 1`

**NOTE:** The script will iterate through all the features (polygons) within the shapefile creating a folder for each specific AOI. For that, you need to specify a column with the 'aoi_name'. In the Git repository you will find an example of a test shapefile.

## Required packages
- Python 3.9
- PDAL
- Geopandas
- wget
- pathlib2

# Installation on specific environment
if you want to create a specific environment, please follow the instructions below

## Python environment with Miniconda

Go to the [Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers) webpage 

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
```

## Installation using mamba
 Install mamba using the next comand line
 ```
 conda install mamba -n base -c conda-forge
 ```
 then, install all packages using the `*.yml` file provided.
 ```
 mamba env create -f py39.yml
 ```
 To active the environment, type `conda activate py39`

## Installation using conda

```
conda env create -f py39.yml
```
To active the environment, type `conda activate py39`