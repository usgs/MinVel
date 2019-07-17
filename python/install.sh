#!/bin/bash

# Adapted O. Boyd from ShakeMap install script, 7/15/19

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    prof=~/.bashrc
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    prof=~/.bash_profile
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
else
    echo "Unsupported environment. Exiting."
    exit
fi

source $prof

# Name of virtual environment
VENV=minvel

developer=0
py_ver=3.6

echo "Using python version $py_ver"

# Is conda installed?
conda --version
if [ $? -ne 0 ]; then
    echo "No conda detected, installing miniconda..."

    command -v curl >/dev/null 2>&1 || { echo >&2 "Script requires curl but it's not installed. Aborting."; exit 1; }

    curl $mini_conda_url -o miniconda.sh;

    # if curl fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to download miniconda installer shell script. Exiting."
        exit 1
    fi
    
    echo "Install directory: $HOME/miniconda"

    bash miniconda.sh -f -b -p $HOME/miniconda

    # if miniconda.sh fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to run miniconda installer shell script. Exiting."
        exit 1
    fi
    
    . $HOME/miniconda/etc/profile.d/conda.sh
else
    . $HOME/miniconda/etc/profile.d/conda.sh
    echo "conda detected, installing $VENV environment..."
fi

echo "Installing packages from conda-forge"

# Choose an environment file based on platform
# only add this line if it does not already exist
grep "/etc/profile.d/conda.sh" $prof
if [ $? -ne 0 ]; then
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $prof
    echo 'export PATH="'$HOME'/miniconda/bin:$PATH"' >> $prof
fi
# Make sure PYTHONPATH is set up correctly
grep "MinVel" $prof
if [ $? -ne 0 ]; then
    echo 'export PYTHONPATH="$PYTHONPATH:'$PWD'"' >> $prof
fi

# Start in conda base environment
conda activate base

# Remove existing minvel environment if it exists
conda remove -y -n $VENV --all

package_list=(
      "python=$py_ver"
      "netCDF4"
      "numpy"
)

# Create a conda virtual environment
echo "Creating the $VENV virtual environment:"
conda create -y -n $VENV -c conda-forge \
      --channel-priority ${package_list[*]}

# Bail out at this point if the conda create command fails.
# Clean up zip files we've downloaded
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment.  Resolve any conflicts, then try again."
    exit
fi

# Activate the new environment
echo "Activating the $VENV virtual environment"
conda activate $VENV

# if conda activate fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to activate ${VENV} conda environment. Exiting."
    exit 1
fi

# upgrade pip, mostly so pip doesn't complain about not being new...
pip install --upgrade pip

# if pip upgrade fails, complain but try to keep going
if [ $? -ne 0 ];then
    echo "Failed to upgrade pip, trying to continue..."
    exit 1
fi

# Tell the user they have to activate this environment
echo ""
echo "Prepare environment with '. "$prof"'"
echo "Type 'conda activate $VENV' to use this new virtual environment."
echo "Test the python installation with 'python tests/test.py'"
echo ""
