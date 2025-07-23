#!/bin/bash

###########################################################
# Download and extract the compiled database for running  #
# Master.                                                 #
#                                                         #
#         Created by Furman Lab at HUJI, 2025.            #
###########################################################

DB_DIR="databases/"
mkdir -p $DB_DIR

############# DATABASE #############
# Download database for Master
echo "Downloading MASTER database. This may take a while..."
MASTERDB="$DB_DIR/masterDB/"
mkdir -p $MASTERDB
curl -o "$DB_DIR/masterDB.tar.gz" -L "https://www.cs.huji.ac.il/~jvarga/master_search.tar.gz"

echo "Extracting MASTER search database..."
tar -xzf "$DB_DIR/masterDB.tar.gz" -C $DB_DIR
find "$(realpath "$DB_DIR/masterDB/")" -name '*pds' > db_list_30nonred # Create list file with full paths

# Download cleaned PDB files for Master
echo "Downloading and extracting cleaned PDBs for MASTER search. This may take a while..."
curl -o "$DB_DIR/master_clean.tar.gz" -L "https://www.cs.huji.ac.il/~jvarga/master_extract.tar.gz"

echo "Extracting cleaned PDBs..."
tar -xzf "$DB_DIR/master_clean.tar.gz" -C "$DB_DIR"

rm "$DB_DIR/master_clean.tar.gz" "$DB_DIR/masterDB.tar.gz"

echo "Downloads are finished!"
