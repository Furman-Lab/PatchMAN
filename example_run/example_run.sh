#!/bin/bash

rec=1N7E_A.pdb
pep=pepseq
native=1N7E_native.pdb
MASTER='path/to/master/bin/'
DB='path/to/DB'
ROSETTA='path/to/rosetta'
PATCHMAN='path/to/patchman'

$PATCHMAN/split_to_motifs.py $rec

ls '???_'$rec` > motif_list

$MASTER/create_pds --type query --pdbList motif_list

ls $DB/*pds > db_list

for a in `ls *pds`; do motif=`echo $a | cut -d '.' -f 1`; $MASTER/master --query $a --targetList db_list --bbRMSD --rmsdCut 1.5 --topN 1000000 --matchOut $motif'_matches'; done

for a in `cat motif_list`; do name=`echo $a | cut -d '.' -f 1`; $PATCHMAN/extract_peps_for_motif.py -m $name'_matches' -p $pep -r $rec --patch $a > $name'_log'; done

ls ???_????_*_*_0001.pdb > input_list

$ROSETTA/main/source/bin/FlexPepDock.linuxgccrelease -in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary -out:file:silent decoys.silent -lowres_preoptimize -flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only -ex1 -ex2aro -use_input_sc -unboundrot $rec -native $native

