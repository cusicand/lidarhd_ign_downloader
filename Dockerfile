FROM debian:bookworm-slim

RUN apt-get update && apt-get -y upgrade && apt-get -y install \
	vim \
	git \ 
	git-core \
	wget \
	unzip \
	cmake \
	gdal-bin \
	libgdal-dev \
	build-essential

WORKDIR /root

# Add the repository
ADD ./* /root/lidarhd_ign_downloader/

# Install conda and dependencies
SHELL ["/bin/bash", "-c"]
RUN mkdir -p ~/miniconda3 && \
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
	bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 && \
	rm ~/miniconda3/miniconda.sh && \
	source ~/miniconda3/bin/activate && \
	conda install -y -n base -c conda-forge mamba pdal python-pdal gdal && \
	cd lidarhd_ign_downloader && \
	mamba env create -f pdal_env.yml

# Set environment variables
ENV HOME="/root"
ENV LIDAR_PATH="${HOME}/lidarhd_ign_downloader"
ENV PATH="${LIDAR_PATH}:${PATH}"            
ENV PYTHONPATH="${PYTHONPATH}:${LIDAR_PATH}"
ENV PATH="${HOME}/miniconda3/bin:$PATH"

# Install PDAL Wrench
RUN cd ~ && \
	git clone https://github.com/PDAL/wrench.git && \
	cd wrench && \
	mkdir build && \
	cd build && \
	cmake .. && \
	make

ENV PDWRENCH="${HOME}/wrench/build/"
ENV PATH="${PDWRENCH}:${PATH}"

# Set the working directory
RUN chmod +x ~/lidarhd_ign_downloader/lidar_downloader.py
WORKDIR /root/lidarhd_ign_downloader
RUN conda init bash \
	&& echo -e "\nconda activate pdal_env" >> ~/.bashrc \
	&& echo "echo \"Lidar HD IGN downloader\"" >>  ~/.bashrc \
	&& echo "echo \"Type lidar_downloader.py -h to access to the help of the tool\"" >>  ~/.bashrc \
	&& conda clean -a

