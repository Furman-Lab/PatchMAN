#!/bin/bash
#SBATCH --time=93:00:00
#SBATCH --mem=2000m

# sort by reweighted_sc
printScoreFile_byHeader.pl score.sc reweighted_sc I_sc rmsBB rmsBB_if description | gawk '{print $2, $3, $4, $5, $6}' > short.sc
sort -nk1 short.sc > sorted.sc

# filter possible disulf bridges
for a in $(head sorted.sc | gawk '{print $NF}'); do score=$(grep $a sorted.sc | gawk '{print $1}' | cut -d '.' -f 1);
 if (( $score < -1000000 ));
  then echo $a >> tmp;
 fi; done;
if [ -f tmp ]; then
all_filter=$(echo `cat tmp` | sed 's/ /|/g');
grep -v -E "$all_filter" sorted.sc > tmp2
mv tmp2 sorted.sc
rm -f tmp
fi

# tale top 1 percent structures
nlines=`wc -l sorted.sc | gawk '{print $1}'`
head -$(($nlines/100)) sorted.sc | gawk '{print $NF}'> top1percent

# run clustering
mkdir clustering
cd clustering
pdb=`ls ../???_????_*_*_0001.pdb | head -1`
len=`grep CA $pdb | wc -l`
peps=$1 # pep seq argument
plen=`echo ${#peps}`
R=$2 #radius arg
actualR=`date | awk '{print sqrt('$plen'/'$len')*'$R'}'`
echo actual radius is "$actualR"
$ROSETTA_BIN/cluster.linuxgccrelease -in:file:silent ../decoys.silent -in:file:silent_struct_type binary -cluster:radius "$actualR" -in:file:fullatom -tags `cat ../top1percent` > clog
