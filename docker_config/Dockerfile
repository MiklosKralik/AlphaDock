# install operating system
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive

# install dependencies
RUN apt-get update -y
RUN apt-get install git build-essential python3-dev libglew-dev \
    libpng-dev libfreetype6-dev libxml2-dev \
    libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev -y
RUN apt-get install pymol -y
RUN apt-get install neovim -y
RUN apt-get install autodock-vina -y
RUN apt-get install openbabel -y
RUN apt-get install apache2-utils  -y
RUN apt-get install openssh-server  -y
RUN apt-get install supervisor -y
RUN apt-get install python3-pip -y
RUN pip install pandas
RUN pip install matplotlib
