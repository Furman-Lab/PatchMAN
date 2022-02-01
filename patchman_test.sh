#!/bin/bash

###############################################
#                Test PatchMAN                #
# This script will run an example job         #
# on the 1ssh structure in the test directory #
# without receptor backbone minimization.     #
# Run from PatchMAN's base directory          #
#                                             #
#     Created by Furman Lab at HUJI, 2022.    #
###############################################

cd test/
../PatchMAN_protocol.sh -m false 1ssh.pdb EGPPPAMPARPT
