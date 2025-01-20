#!/usr/bin/bash
set -e

# Put the URL of the MASTER archive here, after registration here
MASTER_URL=""
PATCHMAN_ENV="$PWD/patchman/" # Can be changed to default installation path: $HOME/micromamba/envs/patchman
MICROMAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest" # adjust this for different OS
MASTER_DIR="$PWD/master"  # Adjust if master is extracted elsewhere

# DO NOT CHANGE ANYTHING BELOW THIS LINE, UNLESS YOU KNOW WHAT YOU ARE DOING

PATCHMAN_ENV_BIN="$PATCHMAN_ENV/bin"

# Download micromamba to manage the environment
curl -fsSL "$MICROMAMBA_URL" | tar -xvj --strip-components=1 -C ./ bin/micromamba
eval "$(./micromamba shell hook -s posix)"

# Install PyRosetta into this directory.
micromamba env create -f environment.yml --prefix "$PATCHMAN_ENV" -y

# Download MASTER and mslib patch
echo "Downloading MASTER and mslib..."
curl -f -o master.tar.gz "${MASTER_URL}"
curl -f -o mslib.tar.gz "https://grigoryanlab.org/msl/msl-static-Linux-x86-64_1.2.2.7.tar.gz"


# Install MASTER into the conda environment
# Step 1: Create a directory for extracting the master archive
mkdir -p "$MASTER_DIR"

# Step 2: Extract the master and mslib archives
echo "Extracting master archive..."
tar -xzf master.tar.gz -C "$MASTER_DIR" --strip-components=1
echo "Extracting mslib archive..."
tar -xzf mslib.tar.gz -C "$PWD/"

# Step 3: Patch the MASTER source code
echo "Patching MASTER source code..."
cp "$MASTER_DIR/src/Match.cpp" "$MASTER_DIR/src/Match_old.cpp" # create backup
dos2unix "$MASTER_DIR/src/Match.cpp"
patch -l "$MASTER_DIR/src/Match.cpp" bin/master.patch

# Step 4: Compile MASTER
echo "Compiling MASTER..."
cd "$MASTER_DIR"
make all
cd ../

# Step 5: Install the compiled executable into the bin directory of patchman_env
echo "Installing MASTER executable..."
mkdir -p "$PATCHMAN_ENV_BIN"
cp "$MASTER_DIR/bin/"* "$PATCHMAN_ENV_BIN/"

# Clean up
echo "Cleaning up..."
cd ..
rm -rf "$MASTER_DIR" master.tar.gz mslib.tar.gz

echo "Installation complete! MASTER is now installed in $PATCHMAN_ENV_BIN"