# IGN LiDAR HD downloader
Python CLI script to download and resample LiDAR HD classified point clouds from IGN-France.

To have some news about data processing, please visit [IGN-France website](https://geoservices.ign.fr/lidarhd).

## Required packages
- python (conda) environment
  - python>=3.10
  - pdal==2.5.6
  - draco=1.5.6 (we should specify draco version because there is an issue on the latest conda python-pdal version)
  - python-pdal
  - untwine
  - geopandas
  - python-wget
  - pathlib2
  - tqdm

## Clone github repository 

In order to install properly, you have to clone `lidar_ign_downloader` repository.
```
git clone https://github.com/cusicand/lidarhd_ign_downloader.git
cd ./lidarhd_ign_downloader
```
## Installation

Most of the libraries required for this script are standard and are often pre-installed in conda python environments. Please follow the instructions below depending on your requirements.

If you already have a conda python environment pre-installed, please follow the instruction in section [Install on pre-existing conda python environment](#install-on-an-existing-python-environment) section. Otherwise, you will need to install a conda python environment to use this script. Instructions are given in the section [Install packages on a new conda python environment](#install-packages-on-a-new-python-environment).

### Install packages on an existing python environment

Run the next command lines in your command-line prompt:

`conda activate <your-env-name>`

`conda install -c conda-forge pdal==2.5.6 draco=1.5.6 python-pdal gdal untwine geopandas python-wget pathlib2`

We encourage the use of `mamba` since this library as is faster than conda. 
If you want to use `mamba`, run the following lines:

`conda activate <your-env-name>`

`mamba install -c conda-forge pdal==2.5.6 draco=1.5.6 python-pdal gdal untwine geopandas python-wget pathlib2`

### Install packages on a new python environment

If you want to create a specific python environment, please follow the instructions below.

#### Python environment with Miniconda

Go to the [Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers) website and download the lastest version of Miniconda. Detailed instructions on how to install conda python environments for your operating system are available on the [Anaconda website](https://docs.anaconda.com/free/miniconda/).

Once conda installed, you can 

#### Create the new environment using `mamba`

First, install `mamba`
```
conda install -n base -c conda-forge mamba 
```
then, install all packages using the `pdal_env.yml` file provided.
```
mamba env create -f pdal_env.yml
```
To active the environment, type `conda activate pdal_env`

#### Create the new environment using `conda`

```
conda env create -f pdal_env.yml
```
To active the environment, type `conda activate pdal_env`

### Make python script executables

If you want to run the script anywhere in your computer from CLI, you need to add the following lines to your `.bashrc` file to have full access to scripts. 
Open your `.bashrc` file using vi `~/.bashrc` or `nano ~/.bashrc` and copy the following lines at the end.

```
export LIDAR_PATH=$HOME/lidarhd_ign_downloader
export PATH=$LIDAR_PATH:$PATH            
export PYTHONPATH=$LIDAR_PATH:$PYTHONPATH
```
Once changes saved, run `chmod +x lidar_downloader.py` into github repository to ensure the execution of the script.

Use `source ~/.bashrc` to reload changes.

## How to use

Iside the command line prompt, type `python lidar_downloader.py -h` to access to the help of the tool.

This small tool needs and Area of Interest (AOI) in shapefile format as mandatory argument.

Then, we can switch between several parameters like:
- 

`some_amazing_name.shp`.
Then, we can switch between several parameters like:
- `-out_data` or `--out_data_path` to select specific path.
- `-tr` or `--dem_resolution` to select the desired resolution of the output DEM.
- `-dtype` or `--file_data_type` to switch between `GTiff` and `VRT` files.

Below is an example using the supplied shapefile:

`lidarhd_downloader.py aoi_double.shp -out_data /home/user/some/path/directory/ -tr 1 -dtype GTiff` or 

`lidarhd_downloader.py aoi_double.shp --out_data_path /home/user/some/path/directory/ --dem_resolution 1 --file_data_type GTiff`

**NOTE:** The script will iterate through all the features (polygons) within the shapefile crating a folder for each specific AOI. For that, you need to specify a column with the `aoi_name`.

## Contact

For any question/bug/issue, please report it on issues section or contact [diego.cusicanqui@univ-grenoble-alpes.fr](mailto:diego.cusicanqui@univ-grenoble-alpes.fr)