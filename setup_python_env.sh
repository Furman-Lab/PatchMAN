#!/bin/bash

###########################################################
#           Setup file for python and PyRosetta           #
#                                                         #
# This script installs pytho packages and PyRosetta. Not  #
# needed if working with containers.                      #
#                                                         #
#           Created by Furman Lab at HUJI, 2022.          #
###########################################################

# Setup virtualenv
virtualenv -p python3 venv_PatchmanProtocol
source venv_PatchmanProtocol/bin/activate

# Install python pacakges and PyRosetta
pip3 install -r requirements.txt
pip3 install containers/pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl
