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
# Register for MASTER at https://grigoryanlab.org/index.php?sec=download&soft=MASTER and insert the URL of MASTER-v1.6 at URL_MASTER in the .env file


source .env
DB_DIR="databases/"
mkdir -p $DB_DIR

############# DATABASE #############
# Download database for Master
echo "Downloading MASTER database. This may take a while..."
MASTERDB="$DB_DIR/masterDB/"
mkdir -p $MASTERDB
rsync -varz arteni.cs.dartmouth.edu::masterDB/ $MASTERDB
find $PWD/$MASTERDB -name '*pds'  > db_list_30nonred # create list file for Master search

# Download cleaned PDB files for Master
echo "Downloading and extracting cleaned PDBs for MASTER search. This may take a while..."
wget -O databases/master_clean.tar.gz https://www.cs.huji.ac.il/~jvarga/master_clean.tar.gz
cd databases/
tar -xzf master_clean.tar.gz

echo "Downloads are finished!"
