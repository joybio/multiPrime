#!/bin/python

"""extract clusters from ch-hit results"""
#caution the chrom column of the position file must contain chr
__date__ = "2022-6-7"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser
import collections
from collections import defaultdict
import re
import os

parser = OptionParser('Usage: %prog ')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: candidate_primers_sets.txt.')
parser.add_option('-o','--out',
                dest='out',
                help='output file: filter candidate primer by primer coverage.')
parser.add_option('-n','--number',
                dest='number',
                default = 10,
                type="int",
                help='min number of primer coverage. Default: 10')
import random
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

(options,args) = parser.parse_args()
data = open(options.input,"r")
out = open(options.out,"w")
number = options.number

for i in data:
	line = i
	i = i.strip().split("\t")
	cluster_info = i[0].split("/")
	cluster_number = int(cluster_info[-1].split("_")[-1].split(".")[0])
	if cluster_number >= number:
		out.write(line)
	else:
		pass
data.close()
out.close()




