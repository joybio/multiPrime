#!/bin/python
"""
Build a binary dict using the columns of input file or fasta title and fasta sequence.
"""

__date__ = "2022-6-1"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import pickle
import re
from collections import defaultdict
from optparse import OptionParser


def argsParse():
	parser = OptionParser(
		'Usage:   %prog -i [input] -n [int] -o [output] -f [txt] options {-v [int or T]} \n'
		'[fasta]: %prog -i [input] -o [output] -f [fa] options {-t [F or T]} \n'
		'It will return a binary dict, which can loaded by pickle.')
	parser.add_option('-i', '--input',
					  dest='input',
					  help='Input file: It can be a txt file or a fasta file.')

	parser.add_option('-n', '--num',
					  dest='num',
					  default="0",
					  type="int",
					  help='Which column should be used as hash key. for example:'
						   'if you want to set column 1 as the key of dict, then -n 0. Default: [0].')

	parser.add_option('-t', '--head',
					  dest='head',
					  default="T",
					  type="str",
					  help='Head information of the fasta, the complete title of fasta means all information after >.'
						   'Default: T. key = Accession number. value = title + fasta.')

	parser.add_option('-f', '--form',
					  dest='form',
					  default="txt",
					  type="str",
					  help='Format of input file. it can be [txt] or [fa]. Default: [txt].')

	parser.add_option('-v', '--value',
					  dest='value',
					  default="T",
					  type="str",
					  help='Value: which column should be used as hash value.'
						   'if value == T(total),use the whole line as the hash value. Default: [T].')

	parser.add_option('-o', '--out',
					  dest='out',
					  help='Out file: binary file, which is a dictionary.')

	(options, args) = parser.parse_args()
	import sys
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


def prepare_pickle(Input, Output):
	pickle_dict = defaultdict(list)
	with open(Input, "r") as In:
		for i in In:
			line = i.strip()
			i = i.strip().split("\t")
			key = i[column]
			if value == "T":
				pickle_dict[key] = line
			else:
				pickle_dict[key].append(i[int(value)])
	with open(Output, "wb") as Out:
		pickle.dump(pickle_dict, Out)


def prepare_fa_pickle(Input, Output):
	pickle_dict = {}
	with open(Input, "r") as In:
		for i in In:
			if i.startswith(">"):
				line = i
				i = i.lstrip(">")
				key_list = []
				if re.search(">", i):
					# >FJ608688.1 ..., partial cds >FJ608691.1 ..., partial cds
					ID_list = i.split(">")
					for j in ID_list:
						j = j.split(" ")
						key_list.append(j[0])
				else:
					i = i.split(" ")
					key_list.append(i[0])
			else:
				if headinfo == "T":
					value = line + i
					for k in key_list:
						pickle_dict[k] = value
				else:
					value = i
					for k in key_list:
						pickle_dict[k] = value
	with open(Output, "wb") as Out:
		pickle.dump(pickle_dict, Out)


if __name__ == "__main__":
	(options, args) = argsParse()
	if options.form == "txt":
		column = int(options.num)
		value = options.value
		prepare_pickle(options.input, options.out)
	elif options.form == "fa":
		headinfo = options.head
		prepare_fa_pickle(options.input, options.out)
