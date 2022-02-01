#!/bin/bash

###########################################################
# Download all files needed for installation of PatchMan, #
# also download and extract the compiled database for     #
# running Master. The indicated versions were used in     #
# the paper.                                              #
#                                                         #
#         Created by Furman Lab at HUJI, 2022.            #
###########################################################


# Get Rosetta and PyRosetta username and password from .env file. Create this file based on sample.env
# For Rosetta, obtain a license from: https://els2.comotion.uw.edu/product/rosetta
# For PyRosetta, obtain a license from: https://els2.comotion.uw.edu/product/pyrosetta


source .env

# Download Rosetta
echo "Downloading Rosetta. This may take a while..."
ROSETTA_LINK="https://www.rosettacommons.org/downloads/academic/2019/wk14/rosetta_src_2019.14.60699_bundle.tgz"
curl -f -o rosetta.tar.gz  -u ${ROSETTA_USER}:${ROSETTA_PASS} ${ROSETTA_LINK} --keepalive-time 2 -H 'Expect:'

# Download PyRosetta
echo "Downloading PyRosetta. This may take a while..."
PYROSETTA_LINK="https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python35.linux.wheel/pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl"
curl -k -f -o pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl -u ${PYROSETTA_USER}:${PYROSETTA_PASS} ${PYROSETTA_LINK} --keepalive-time 2 -H 'Expect:'

# Download MASTER and mslib patch
echo "Downloading MASTER and mslib..."
curl -f -o master.tar.gz https://grigoryanlab.org/master/master-v1.6.tar.gz
curl -f -o mslib.tar.gz https://grigoryanlab.org/msl/msl-static-Linux-x86-64_1.2.2.7.tar.gz

# Download database for Master
echo "Downloading and extracting database for Master search. This may take a while..."
wget masterDB.tar.gz https://www.cs.huji.ac.il/~jvarga/masterDB.tar.gz
tar -xzf masterDB.tar.gz
find $PWD/masterDB/ -name '*pds'  > db_list_30nonred # create list file for Master search

# Download OpenMPI that matches with host distribution. If you need an other one, change OMPI_VERSION.
echo " Downloading OpenMPI..."
OMPI_VERSION=$(mpirun --version | head -1 | cut -d ' ' -f4) # change this to match your mpi version!
OMPI_URL="https://download.open-mpi.org/release/open-mpi/v2.1/openmpi-$OMPI_VERSION.tar.bz2"
curl -Lo openmpi.tar.bz2 $OMPI_URL --keepalive-time 2 -H 'Expect:'
