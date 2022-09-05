#!/bin/pyhthon
"""genome format"""
#caution the chrom column of the position file must contain chr
__date__ = "2022-6-1"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input.fa] -o [output.fa]')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: Check')

parser.add_option('-o','--out',
                dest='out',
                help='Out file: fasta')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
import os 
os.system("cat {} | grep {} > {}".format(options.input,,))
os.system("bedtools getfasta -fi {} -bed {} -fo {}".format())



