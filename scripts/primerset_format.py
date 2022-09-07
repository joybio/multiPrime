#!/bin/python
"""format"""
# equal to
# awk -F '/' '{print $NF}' final_maxprimers_set.xls | sed 's/.candidate.primers.txt//g' | awk '{print ">"$1"_F\n"$2"\n>"$1"_R\n"$3}'
# caution the chrom column of the position file must contain chr
__date__ = "2022-7-28"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

from optparse import OptionParser


def argsParse():
	parser = OptionParser('Usage: %prog -i [input.xls] -o [output.fa] \n')

	parser.add_option('-i', '--input',
					  dest='input',
					  help='Input file: final_maxprimers_set.xls')

	parser.add_option('-o', '--out',
					  dest='out',
					  help='Out file: final_maxprimers_set.fa')

	(options, args) = parser.parse_args()
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


def primerset_format(Input, Output):
	with open(Input, "r") as In:
		with open(Output,"w") as out:
			for i in In:
				if i.startswith("#"):
					pass
				else:
					i = i.strip().split("/")
					info = i[-1].replace(".candidate.primers.txt", "").split("\t")
					out.write(">" + info[0] + "_F\n" + info[2] + "\n" + \
							  ">" + info[0] + "_R\n" + info[3] + "\n")


if __name__ == "__main__":
	(options, args) = argsParse()
	primerset_format(options.input, options.out)
