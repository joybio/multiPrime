#!/bin/bash
if [ $# != 2 ];then
echo -e "primer_filter.sh <new primer.fa> <core primer set.fa>"
exit
fi
new=$1
core=$2

python /share/data3/yangjunbo/git_storage/multiPrime/scripts/Primer_set_update.py -c $core -n $new -r /share/data3/yangjunbo/database/T2T_gtf/T2T_bowtie2.len300 -p 10 -f D -o ${new%.fa}

awk '{if($1!~/Cluster/){print $8"\t"$1}else if($8!~/Cluster/){print $1"\t"$8}}' ${new%.fa}.dimer | sort > ${new%.fa}.dimer.filter.sort

sed 'N;s/\n/\t/g' $new | sort | uniq > ${new%.fa}.format.sort
join -a1 ${new%.fa}.format.sort ${new%.fa}.dimer.filter.sort | awk -F ' ' '$3!~/>/' > ${new%.fa}.update.fa


