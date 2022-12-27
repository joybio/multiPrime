#!/bin/bash
if [ $# != 2 ];then
echo -e "selct.sh <query fa / primer set> <dimer out>
# select primer by dimerzation"
exit
fi
fa=$1 # candidate priemr file
pool=$2 # primer pool file

cat $fa $pool > merge.primer.pool.fa
python /share/data3/yangjunbo/git_storage/multiPrime/scripts/finDimer.py -i merge.primer.pool.fa -n 20 -o merge.primer.pool.fa.dimer
sed 'N;s/\n/\t/g' $fa | sort > candidate.primer.sort.txt
sort merge.primer.pool.fa.dimer > merge.primer.pool.fa.dimer.sort
join -a1 -a2 -o 1.1 1.2 2.1 candidate.primer.sort.txt merge.primer.pool.fa.dimer.sort | awk 'NF==2' | sort | uniq  > candidate.primer.fa

