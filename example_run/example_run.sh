#!/bin/bash

rec=1OOT_A.pdb
pep=1SSH_B.fasta
native=1OOT_native.pdb
MASTER='/vol/ek/Home/alisa/tools/master/master-v1.6/bin'
DB='path/to/DB/'
ROSETTA='/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle'
PATCHMAN='/vol/ek/Home/alisa/scripts/patchman'

#Split surface to structural patches
$PATCHMAN/split_to_motifs.py $rec

# Create a list of motifs
ls ???'_'$rec > motif_list

# Create pds files for running MASTER
$MASTER/createPDS --type query --pdbList motif_list

# Create a list of database structures for the template search
ls $DB/??/*pds > db_list

# Run MASTER for all motifs
for a in `ls *pds`; do motif=`echo $a | cut -d '.' -f 1`; $MASTER/master --query $a --targetList db_list --bbRMSD --rmsdCut 1.5 --topN 1000000 --matchOut $motif'_matches'; done

# Prepack the receptor structure for further FlexPepDock refinement
$ROSETTA/main/source/bin/FlexPepDocking.linuxgccrelease -s $native -out:pdb -scorefile ppk.score.sc -nstruct 1 -flexpep_prepack -ex1 -ex2aro -use_input_sc -unboundrot $rec

# Suggesting that the receptor chain is chain A:

rec_name=`echo $rec | cut -d '.' -f 1`
grep ' A ' *native_0001.pdb > $rec_name'.ppk.pdb'

# Extract templates and thread the peptide sequence
for a in `cat motif_list`; do name=`echo $a | cut -d '.' -f 1`; $PATCHMAN/extract_peps_for_motif.py -m $name'_matches' -p $pep -r $rec_name'.ppk.pdb' --patch $a > $name'_log'; done

# Create a list of input structures for refinement
ls ???_????_*_*_0001.pdb > input_list

# Run FlexPepDock refinement
$ROSETTA/main/source/bin/FlexPepDocking.linuxgccrelease -in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary -out:file:silent decoys.silent -min_receptor_bb -lowres_preoptimize -flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only -ex1 -ex2aro -use_input_sc -unboundrot $rec -native $native

