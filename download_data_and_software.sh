#!/bin/bash

###########################################################
# Download all files needed for installation of PatchMan, #
# also download and extract the compiled database for     #
# running Master. The indicated versions were used in     #
# the paper.                                              #
#                                                         #
#         Created by Furman Lab at HUJI, 2024.            #
###########################################################


# Get Rosetta and PyRosetta username and password from .env file. Create this file based on sample.env
# For Rosetta, obtain a license from: https://els2.comotion.uw.edu/product/rosetta
# For PyRosetta, obtain a license from: https://els2.comotion.uw.edu/product/pyrosetta
# Register for MASTER at https://grigoryanlab.org/index.php?sec=download&soft=MASTER and insert the URL of MASTER-v1.6 at URL_MASTER in the .env file


source .env
CONTAINER_DIR="containers/"
DB_DIR="databases/"
mkdir -p $CONTAINER_DIR $DB_DIR

############# SOFTWARE #############
# Download Rosetta
echo "Downloading Rosetta. This may take a while..."
ROSETTA_LINK="https://www.rosettacommons.org/downloads/academic/2019/wk14/rosetta_src_2019.14.60699_bundle.tgz"
curl -f -o $CONTAINER_DIR/rosetta.tar.gz  -u ${ROSETTA_USER}:${ROSETTA_PASS} ${ROSETTA_LINK} --keepalive-time 2 -H 'Expect:'

# Download PyRosetta
echo "Downloading PyRosetta. This may take a while..."
PYROSETTA_LINK="https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python35.linux.wheel/pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl"
curl -k -f -o $CONTAINER_DIR/pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl -u ${PYROSETTA_USER}:${PYROSETTA_PASS} ${PYROSETTA_LINK} --keepalive-time 2 -H 'Expect:'

# Download MASTER and mslib patch
echo "Downloading MASTER and mslib..."
curl -f -o $CONTAINER_DIR/master.tar.gz "${MASTER_URL}"
curl -f -o $CONTAINER_DIR/mslib.tar.gz "https://grigoryanlab.org/msl/msl-static-Linux-x86-64_1.2.2.7.tar.gz"

# Download OpenMPI that matches with host distribution. If you need an other one, change OMPI_VERSION.
echo " Downloading OpenMPI..."
OMPI_VERSION=$(mpirun --version | head -1 | cut -d ' ' -f4) # change this to match your mpi version!
OMPI_SHORT_VERSION=$(echo $OMPI_VERSION | cut -d '.' -f1-2)
OMPI_URL="https://download.open-mpi.org/release/open-mpi/v$OMPI_SHORT_VERSION/openmpi-$OMPI_VERSION.tar.bz2"
curl -Lo $CONTAINER_DIR/openmpi.tar.bz2 $OMPI_URL --keepalive-time 2 -H 'Expect:'



############# DATABASE #############
# Download database for Master
echo "Downloading MASTER database. This may take a while..."
MASTERDB="$DB_DIR/masterDB/"
mkdir -p $MASTERDB
wget https://zenodo.org/records/13118411/files/master_search.tar.gz?download=1 -o $MASTERDB.tar.gz
find $PWD/$MASTERDB -name '*pds'  > db_list_30nonred # create list file for Master search

# Download cleaned PDB files for Master
echo "Downloading and extracting cleaned PDBs for MASTER search. This may take a while..."
wget -O databases/master_clean.tar.gz https://zenodo.org/records/13118411/files/master_extract.tar.gz?download=1
cd databases/
tar -xzf master_clean.tar.gz

echo "Downloads are finished!"
