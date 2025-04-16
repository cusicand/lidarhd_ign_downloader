# IGN LiDAR HD downloader

[![DOI](https://zenodo.org/badge/706232299.svg)](https://zenodo.org/doi/10.5281/zenodo.10697626)

**Author**: Diego Cusicanqui

**Contributors**: Jean Baptiste Barré

**contact**: [diego.cusicanqui@univ-grenoble-alpes.fr](mailto:diego.cusicanqui@univ-grenoble-alpes.fr)

Python scripts to quickly download and resample LiDAR HD classified point clouds from IGN-France.  
This Python script can be executed on the command-line interface (CLI) as well as within a Jupyter Notebook. Refer to the [How to use section](#how-to-use) for more details.

Visit the [IGN-France website](https://geoservices.ign.fr/lidarhd) for updates on data processing.

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

Most of the libraries required for this script are standard and are often pre-installed in conda Python environments. Please follow the instructions below based on your requirements.

If you already have a conda Python environment pre-installed, please follow the instructions in the section [Install on an existing conda Python environment](#install-on-an-existing-python-environment). Otherwise, you will need to install a conda Python environment to use this script. Instructions are provided in the section [Install packages on a new conda Python environment](#install-packages-on-a-new-python-environment).

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

After version 2.0 of `lidar_downloader`, we introduced the possibility to generate a point density map to quantitatively estimate the number of points at a given resolution. This task is based on `pdal_wrench` and requires a separate installation. A detailed description can be found on the [wrench GitHub](https://github.com/PDAL/wrench). Further investigations are ongoing to better integrate `pdal_wrench` within `lidar_downloader`.

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

Inside the command line prompt, type `python lidar_downloader.py -h` to access the help for the tool.

This tool requires an Area of Interest (AOI) in *.shp or *.gpkg format as a mandatory argument. After version 3.0, several AOIs can be specified. Extracted DEM s=will take the name of the AOI file by default.  Additionally, several parameters can be specified, such as:

- `-out_data` or `--out_data_path`: Specifies the output data directory. If not provided, data will be stored in the base path of `lidar_downloader.py` by default.
- `-tr` or `--dem_resolution`: Sets the desired resolution of the output DEM.
- `-compute_elev` or `--compute_elevation`: Determines how elevation is computed. Options are `{mean, min, max, median}`. The default value is `mean`. This parameter is still experimental.
- `-dtype` or `--file_data_type`: Switches between `gtif` and `vrt` file formats.
- `-force_database` or `--force_redownload_database`: After version 3.0, we manage to download database from WFS service. Use this parameter re-download LiDAR HD database. Previous file will be deleted. Default value: `False`.
- `-rm_tiles` or `--remove_tiles`: Removes individual downloaded tiles after processing. Default value: `True`.
- `-pdensity` or `--point_density_map`: Generates a point density map for the given resolution. **Requires `pdal_wrench` to be installed.** Default value: `True`. Refer to the [PDAL_WRENCH installation (optional)](#pdal_wrench-installation-optional) section for more details.
- `-cpu_w` or `--cpu_workload`: Specifies the CPU usage ratio for processing data using a multi-threaded strategy. This value represents the percentage of CPU resources allocated for processing. **Warning:** It is recommended to reserve some CPUs for the OS and other processes (at least 4 CPUs). Default value: `0.6`. Maximum value: `0.8`.

#### Recommended use
Below is an example using the supplied shapefile:
```bash
lidar_downloader.py /path/to/aoi_example.gpkg -out_data /home/user/some/path/directory/ -tr 1 -compute_elev mean -dtype gtif -rm_tiles -pdensity -cpu_w 0.6
```
or with several AOI files:
```bash
lidar_downloader.py /path/to/aoi_example1.gpkg /path/to/aoi_example2.gpkg /path/to/aoi_example2.gpkg --out_data_path /home/user/some/path/directory/ --dem_resolution 1 --compute_elevation mean --file_data_type gtif --remove_tiles --point_density_map --cpu_workload 0.6
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
> After version 3.0, the script will iterate through all the *.gpkg files provided in the command line [see recommended use)](#recommended-use). It will create a folder for each specific AOI based on their filename by default. AOI file can contain any projection (EPSG:4326).

# Contact and citation ![DOI](https://zenodo.org/badge/706232299.svg)
For any question/bug/issue regarding this tool, please report it on issues section or contact [diego.cusicanqui@univ-grenoble-alpes.fr](mailto:diego.cusicanqui@univ-grenoble-alpes.fr).

Regarding the LiDAR HD data, please visit their [webpage](https://geoservices.ign.fr/lidarhd) to know more about how to cite their data.

> [!IMPORTANT]   
> If you use this tool, please cite using the following [DOI](https://zenodo.org/doi/10.5281/zenodo.10697626). This will allow some recognition of the time invested and open access to this tool. 
![DOI](https://zenodo.org/badge/706232299.svg)