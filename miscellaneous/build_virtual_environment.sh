#!/bin/bash

VENV_HOME=$PWD/venv3
# Prepare local environment:
python3 -m venv $VENV_HOME

# Activate python3 environment:
source $VENV_HOME/bin/activate

# Update pip and install modules:
pip install --upgrade pip
pip install astropy
# Could need tkinter on box (ex: dnf install python3-tkinter)
pip install matplotlib seaborn
pip install defusedxml
pip install numpy scipy
pip freeze

# more libs
# pip install numpy ctools gammalib gammapy setuptools

# requirements:
#   gcc
#   gcc-c++
#   cfitsio
#   cfitsio-devel
#   doxygen  (optional)
#   readline-devel
function install_ctools_lib() {
	local DOWNLOADS=$PWD/downloads
	mkdir $DOWNLOADS
	cd $DOWNLOADS
	wget http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.6.2.tar.gz
	wget http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.6.2.tar.gz
	tar xzf ctools-1.6.2.tar.gz
	tar xzf gammalib-1.6.2.tar.gz
	cd $DOWNLOADS/gammalib-1.6.2
	./configure --prefix=$VENV_HOME
	make install
	cd $DOWNLOADS/ctools-1.6.2
	./configure --prefix=$VENV_HOME
	make install
}
# install ctools:
# install_ctools_lib

source ctools_activation.sh $VENV_HOME
# to get the ctools activation:
# export GAMMALIB=$VENV_HOME
# export CTOOLS=$VENV_HOME
# source $GAMMALIB/bin/gammalib-init.sh
# source $CTOOLS/bin/ctools-init.sh
