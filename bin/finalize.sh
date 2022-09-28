#!/bin/bash
FINAL_RESULTS='results'

cd clustering;
x=`wc -l ../top1percent | awk '{print $1}'`;
tail -$((x+5)) clog | head -${x} | grep _0001 | awk  '{for(i=1;i<=NF;i++){if ($i ~ /_0001/){print $i,$(i+1),$(i+2)}}}' | sed '/^$/d' > cluster_list;
for i in `awk '{print $1}' cluster_list`;
  do sed 1d ../score.sc | head -1 > tmp1;
    grep $i ../score.sc >> tmp1;
    /vol/ek/share/labscripts/printScoreFile_byHeader.pl tmp1 I_sc reweighted_sc rmsBB_if | tail -1 | awk '{print $2,$3,$4}' >>pdb_list_sc;
  done;
paste cluster_list pdb_list_sc >cluster_list_sc;
echo "Decoy_ID Cluster_no Member_ID I_sc reweighted_sc rmsBB_if" >cluster_list_I_sc_sorted;
echo "Decoy_ID Cluster_no Member_ID I_sc reweighted_sc rmsBB_if" >cluster_list_reweighted_sc_sorted;
sort -nk 4 cluster_list_sc | sort -u -k2,2 | sort -nk 4 | head -10 >>cluster_list_I_sc_sorted;
sort -nk 5 cluster_list_sc | sort -u -k2,2 | sort -nk 5 | head -10 >>cluster_list_reweighted_sc_sorted;
cd ../

mkdir ${FINAL_RESULTS}
cp clustering/cluster_list_reweighted_sc_sorted "${FINAL_RESULTS}/final_scores"
cp "${FINAL_RESULTS}/final_scores" results_info
for a in $(head -11 clustering/cluster_list_reweighted_sc_sorted | tail | gawk '{print "c."$2"."$3".pdb"}');
  do cp clustering/"$a" ${FINAL_RESULTS};
done;

tar -cvzf results.tar.gz $FINAL_RESULTS
