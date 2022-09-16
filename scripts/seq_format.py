#!/bin/python
"""genome format"""
#caution the chrom column of the position file must contain chr
__date__ = "2022-6-1"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser
def argsParse():
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
	elif options.input is None:
		parser.print_help()
		print("Input file must be specified !!!")
		sys.exit(1)
	elif options.out is None:
		parser.print_help()
		print("No output file provided !!!")
		sys.exit(1)
	return parser.parse_args()

def seq_format(Input,Output):
	with open(Input,"r") as In:
		head = In.readline()
		with open(Output,"w") as Out:
			Out.write(head)
			for i in In:
				if i.startswith(">"):
					Out.write("\n" + i)
				else:
					i = i.strip()
					Out.write(i)
			Out.write("\n")

if __name__ == "__main__":
	(options, args) = argsParse()
	seq_format(options.input,options.out)


