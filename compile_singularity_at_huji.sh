rundeb10 -root -user -bind $PWD -run "cd $PWD/containers; singularity build rosetta.sif rosetta.def; singularity build master.sif master.def; singularity build pyrosetta.def pyrosetta.sif"
