#!/bin/python
"""genome format"""
#caution the chrom column of the position file must contain chr

__date__ = "2022-10-8"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

"""
The MIT License (MIT)

Copyright (c) 2022 Junbo Yang <yang_junbo_hi@126.com> <1806389316@pku.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import optparse
from optparse import OptionParser
from collections import defaultdict
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
	parser.add_option('-c','--c',
                        dest='complete',
			default="T",
			type='str',
                        help='Only complete CDS or genome is used, defalut: T. use other word (F), if you dont want this param.')
	parser.add_option('-l','--length',
			dest='length',
			default=300,
			type="int",
			help='Filter fasta by length. Default: 250.')


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

TRANS = str.maketrans("U", "T")
import re

def Trans(seq):
	return seq.translate(TRANS)
seq_dict = defaultdict(str)
seq_lenght_dict = defaultdict(int)
def seq_format(Input):
	complete_number = 0
	with open(Input,"r") as In:
		for i in In:
			if i.startswith(">"):
				key = i
				if re.search("complete", key):
					complete_number += 1
			elif i == "--\n":
				pass
			else:
				i = i.strip()
				seq_dict[key] += i.replace("N","")
				seq_lenght_dict[key] += len(i)
	In.close()
	return seq_dict, seq_lenght_dict, complete_number

if __name__ == "__main__":
	(options, args) = argsParse()
	seq, length, c_number = seq_format(options.input)
	with open(options.out,"w") as Out:
		# if ID dont have string complete, then use all.
		if c_number == 0:
			for i in length.keys():
				if length[i] < options.length:
					pass
				else:
					Out.write(i + seq[i] + "\n")
		else:
			if options.complete == "T":
				for i in length.keys():
					if re.search('complete', i):
						if length[i] < options.length:
							pass
						else:
							Out.write(i + seq[i] + "\n")
			else:
				for i in length.keys():
					if length[i] < options.length:
						pass
					else:
						Out.write(i + seq[i] + "\n")
	Out.close()

