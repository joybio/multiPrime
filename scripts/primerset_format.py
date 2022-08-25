#!/bin/python
"""format"""
# equal to
# awk -F '/' '{print $NF}' final_maxprimers_set.xls | sed 's/.candidate.primers.txt//g' | awk '{print ">"$1"_F\n"$2"\n>"$1"_R\n"$3}'
# caution the chrom column of the position file must contain chr
__date__ = "2022-7-28"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input.xls] -o [output.fa] \n')

parser.add_option('-i','--input',
		dest='input',
		help='Input file: final_maxprimers_set.xls')

parser.add_option('-o','--out',
		dest='out',
		help='Out file: final_maxprimers_set.fa')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
out = open(options.out,"w")
with open(options.input,'r') as f:
	for i in f:
		if i.startswith("#"):
			pass
		else:
			i = i.strip().split("/")
			info = i[-1].replace(".candidate.primers.txt","").split("\t")
			print(info)
			out.write(">" + info[0] + "_F\n" + \
				info[2] + "\n" + \
				">" + info[0] + "_R\n" + \
				info[3] + "\n")
out.close()

