# IGN LiDAR HD downloader

[![DOI](https://zenodo.org/badge/706232299.svg)](https://zenodo.org/doi/10.5281/zenodo.10697626)

**Author**: Diego Cusicanqui

**Contributors**: Jean Baptiste Barré

**contact**: [diego.cusicanqui@univ-grenoble-alpes.fr](mailto:diego.cusicanqui@univ-grenoble-alpes.fr)

Python scripts to quick download and resample LiDAR HD classified point clouds from IGN-France.
This python script can be executed on command line interface (CLI) as well as within a Jupyter Notebook. Look at the [How to use section](#how-to-use) for more details.

Refer to [IGN-France website](https://geoservices.ign.fr/lidarhd) for news about data processing.

## Contents

1. [Required packages](#required-packages)
2. [Clone repository](#clone-github-repository)
3. [Installation](#installation)
4. [Pdal_wrench installation (optional but very useful)](#pdal_wrench-installation-optional)
5. [How to use](#how-to-use)
6. [Contact and citation](#contact-and-citation)

## Required packages
- python (conda) environment
  - python>=3.10
  - pdal==2.5.6
  - draco=1.5.6 (we should specify draco version because there is an issue on the latest conda python-pdal version)
  - python-pdal
  - gdal
  - untwine
  - geopandas
  - python-wget
  - pathlib2
  - tqdm
  - ipykernel
  - pyogrio
  - pdal_wrench (optional if point density map is desired)

## Clone github repository 

In order to install properly, you have to clone `lidar_ign_downloader` repository.
```bash
git clone https://github.com/cusicand/lidarhd_ign_downloader.git
cd ./lidarhd_ign_downloader
```
## Installation

Most of the libraries required for this script are standard and are often pre-installed in conda python environments. Please follow the instructions below depending on your requirements.

If you already have a conda python environment pre-installed, please follow the instruction in section [Install on pre-existing conda python environment](#install-on-an-existing-python-environment) section. Otherwise, you will need to install a conda python environment to use this script. Instructions are given in the section [Install packages on a new conda python environment](#install-packages-on-a-new-python-environment).

### Install packages on an existing python environment

Run the next command lines in your command-line prompt:

```bash
conda activate <your-env-name>
```
or 
```bash
conda install -c conda-forge pdal==2.5.6 draco=1.5.6 python-pdal gdal untwine geopandas python-wget pathlib2 tqdm ipykernel pyogrio
```

We encourage the use of `mamba` since this library as is faster than conda. 
If you want to use `mamba`, run the following lines:

```bash
conda activate <your-env-name>
```

```bash
mamba install -c conda-forge pdal==2.5.6 draco=1.5.6 python-pdal gdal untwine geopandas python-wget pathlib2 tqdm
```

### Install packages on a new python environment

If you want to create a specific python environment, please follow the instructions below.

#### Python environment with Miniconda

Go to the [Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers) website and download the lastest version of Miniconda. Detailed instructions on how to install conda python environments for your operating system are available on the [Anaconda website](https://docs.anaconda.com/free/miniconda/).

Once conda installed, you can 

#### Create the new environment using `mamba`

First, install `mamba`
```bash
conda install -n base -c conda-forge mamba 
```
then, install all packages using the `pdal_env.yml` file provided.
```bash
mamba env create -f pdal_env.yml
```
To active the environment, type `conda activate pdal_env`

#### Create the new environment using `conda`

```bash
conda env create -f pdal_env.yml
```
To active the environment, type:

```bash
conda activate pdal_env
```

### Make python script executables

If you want to run the script anywhere in your computer from CLI, you need to add the following lines to your `.bashrc` file to have full access to scripts. 
Open your `.bashrc` file using vi `~/.bashrc` or `nano ~/.bashrc` and copy the following lines at the end.

```bash
export LIDAR_PATH=$HOME/lidarhd_ign_downloader
export PATH=$LIDAR_PATH:$PATH            
export PYTHONPATH=$LIDAR_PATH:$PYTHONPATH
```
> [!NOTE]   
> If your installation directory is different than `$HOME`, Replace `$HOME` by the full directory path.

Once changes saved, run `chmod +x lidar_downloader.py` into github repository to ensure the execution of the script.

Use `source ~/.bashrc` to reload changes.

## PDAL_WRENCH installation (optional but very useful)

After version 2.0 of `lidar_downloader` we introduce the possibility to generate point density map to quantitatively estimate the number of points at a given resolution. This tasks is based on `pdal_wrench` and requires individual installation. A detailed description could be found [wrench GitHub](https://github.com/PDAL/wrench). Further investigation are ongoing to better integrate pdal_wrench within `lidar_downloader`.

> [!NOTE]
> `Pdal_wrench density` tool has been tested in high mountain environments where vegetation is not an issue. Future test should be conducted in environments with denser vegetation.
> 
### Make pdal_wrench scripts executables

Once `pdal_wrench` installed, it is better you can run it from everywhere in the computer. To do so, add the following lines to your `.bashrc`.
```bash
export PDWRENCH=$HOME/wrench/build/
export PATH=$PDWRENCH:$PATH
```
> [!NOTE]
> If your installation directory is different than `$HOME`, Replace `$HOME` by the full directory path. 

## HOW TO USE

### Command line interface (CLI)

Inside the command line prompt, type `python lidar_downloader.py -h` to access to the help of the tool.

This small tool needs and Area of Interest (AOI) in *.shp or *.gpkg format as mandatory argument. Then, we can switch between several parameters like:

Then, we can switch between several parameters like:
- `-out_data` or `--out_data_path` Out data-path directory. If not specified, data will be stored in lidar_downloader.py base-path by default.
- `-tr` or `--dem_resolution` to select the desired resolution of the output DEM.
- `-compute_elev` or `--compute_elevation` to select the way of how compute elevation. options are {mean,min,max,median}. Default value : mean. This parameter is still experimental.
- `-dtype` or `--file_data_type` to switch between `GTiff` and `VRT` files.
- `-force_database` or  `--force_redownload_database` to force re-download of IGN database. Default value: False.
- `-rm_tiles` or `--remove_tiles` to remove individual downloaded tiles after processing. Default value: True.
- `-pdensity` or `--point_density_map` to generates point density map for given resolution. **Requires pdal_wrench installed**. Default value: True. Refer to [PDAL_WRENCH installation (optional)](#pdal_wrench-installation-optional).
- `-cpu_w` or `--cpu_workload` ratio allow to process the data using a multi-thread strategy. This value represents the percentage of CPU used for processing. **Warning:** Its recommended to keep some CPUs for the OS and other processes (at least 4 CPUs). Default value: 0.6. Max value = 1.0.

#### Recommended use
Below is an example using the supplied shapefile:
```bash
lidarhd_downloader.py aoi_example.gpkg -out_data /home/user/some/path/directory/ -tr 1 -compute_elev mean -dtype gtif -rm_tiles -pdensity -cpu_w 0.6
```
or 
```bash
lidarhd_downloader.py aoi_double.shp --out_data_path /home/user/some/path/directory/ --dem_resolution 1 --compute_elevation mean --file_data_type gtif --remove_tiles --point_density_map --cpu_workload 0.6
```

### Jupyter notebook interface

Inside a Jupyter-Notebook, you can also run the `lidar_downloader`. However, the set up of parameters are is slightly different. You can see and use the [example](./lidar_downloader_jupyter.ipynb) provided within the project.

```python
import lidar_downloader
```

```python
args_list = ['/path/to/aoi/aoi_example.gpkg',
             '--out_data_path', '/some/path/directory/lidar_ign_test/',
             '--dem_resolution', '1.0',
             '--compute_elevation', 'mean',
             '--file_data_type', 'gtif',
             '--remove_tiles',
             '--point_density_map',
             '--cpu_workload', '0.6']
```

```python
parser = lidar_downloader.getparser()
args = parser.parse_args(args_list)
```

```python
lidar_downloader.main(args)
```

> [!NOTE]   
> Whatever the case, the script will iterate through all the features (polygons) within the shapefile or geopackage file. It will create a folder for each specific AOI based on the column with the `aoi_name`. If you used your own shapefile, make sure to have one column called `aoi_name`. Otherwise, you can edit the provided file.

# Contact and citation ![DOI](https://zenodo.org/badge/706232299.svg)
For any question/bug/issue regarding this tool, please report it on issues section or contact [diego.cusicanqui@univ-grenoble-alpes.fr](mailto:diego.cusicanqui@univ-grenoble-alpes.fr).

Regarding the LiDAR HD data, please visit their [webpage](https://geoservices.ign.fr/lidarhd) to know more about how to cite their data.

> [!IMPORTANT]   
> If you use this tool, please cite using the following [DOI](https://zenodo.org/doi/10.5281/zenodo.10697626). This allows to give some recognition for the time invested and open access to this tool. 
![DOI](https://zenodo.org/badge/706232299.svg)