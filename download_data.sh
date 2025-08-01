#!/bin/bash
set -eu
set -o pipefail

hf_url() {
	filename="$1"
  echo "https://huggingface.co/datasets/furman-lab/patchman_datasets/resolve/main/${filename}?download=true"
}

###########################################################
# Download and extract the compiled database for running  #
# Master.                                                 #
#                                                         #
#         Created by Furman Lab at HUJI, 2025.            #
###########################################################

DB_DIR="databases/"
mkdir -p "$DB_DIR"

############# DATABASE #############
# Download database for Master
echo "Downloading MASTER database. This may take a while..."
MASTERDB="$DB_DIR/masterDB/"
mkdir -p $MASTERDB
wget --continue --timeout=0 --tries=10 --retry-connrefused --waitretry=5  \
--user-agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36" -O "$DB_DIR/masterDB.tar.gz" $(hf_url master_search.tar.gz)

echo "Extracting MASTER search database..."
tar -xzf "$DB_DIR/masterDB.tar.gz" -C $DB_DIR
find "$(realpath "$DB_DIR/masterDB/")" -name '*pds' > db_list_30nonred # Create list file with full paths
rm "$DB_DIR/masterDB.tar.gz" -f

# Download cleaned PDB files for Master
echo "Downloading and extracting cleaned PDBs for MASTER search. This may take a while..."
wget --continue --timeout=0 --tries=10 --retry-connrefused --waitretry=5  \
--user-agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36" -O "$DB_DIR/master_clean.tar.gz" $(hf_url master_extract.tar.gz)

echo "Extracting cleaned PDBs..."
tar -xzf "$DB_DIR/master_clean.tar.gz" -C "$DB_DIR"

rm "$DB_DIR/master_clean.tar.gz"

echo "Downloads are finished!"
