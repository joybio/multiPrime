#!/bin/python
"""genome format"""
#caution the chrom column of the position file must contain chr
__date__ = "2022-6-1"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input.fa] -o [output.fa] \n \
	Before: \n \
	> ID1 \n \
	AATTCTTTCTTATC \n \
	CTTCATCTATCATC \n \
	GCAGTCTACGTACG \n \
	...\n \n \
	After: \n \
	> ID1\n \
	AATTCTTTCTTATCCTTCATCTATCATCGCAGTCTACGTACG...')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: raw.fa')
parser.add_option('-o','--out',
                dest='out',
                help='Out file: format.fa')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

data = open(options.input,"r")
out = open(options.out,"w")
head = data.readline()
out.write(head)
for i in data:
	if i.startswith(">"):
		out.write("\n"+i)
	else:
		i = i.strip()
		out.write(i)
out.write("\n")
data.close()
out.close()

