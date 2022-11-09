#!/bin/bash
# select candidate primers for core primerset from new primer set. dependency: findimer or Primer_set_update.py
if [ $# != 3 ];then
echo -e "primer_filter.sh <Path to Primer_set_update.py> <new primer.fa> <core primer set.fa>"
exit
fi
script=$1
new=$2
core=$3

python $script -c $core -n $new -p 10 -f D -o ${new%.fa}

awk '{if($1!~/Cluster/){print $8"\t"$1}else if($8!~/Cluster/){print $1"\t"$8}}' ${new%.fa}.dimer | sort > ${new%.fa}.dimer.filter.sort

sed 'N;s/\n/\t/g' $new | sort | uniq > ${new%.fa}.format.sort
join -a1 ${new%.fa}.format.sort ${new%.fa}.dimer.filter.sort | awk -F ' ' '$3!~/>/' > ${new%.fa}.update.fa


