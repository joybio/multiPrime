#!/bin/python
"""mafft --auto in > out"""
#caution the chrom column of the position file must contain chr

__date__ = "2023-2-24"
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
import os
import optparse
import time
import re
from optparse import OptionParser
from collections import defaultdict
def argsParse():
	parser = OptionParser('Usage: %prog -i [input.fa] -o [output.fa] -s [script]'
				'options: {-l 18 -d 4}')
	parser.add_option('-i','--input',
			dest='input',
			help='Input file')
	parser.add_option('-o','--out',
			dest='out',
			help='Out file')
	parser.add_option('-s','--script',
                        dest='script',
                        help='dir of script')
	parser.add_option('-l','--length',
                        dest='length',
			default='18',
			type="int",
                        help='primer length')
	parser.add_option('-d','--deg',
                        dest='deg',
			default='4',
			type="int",
                        help='primer degeneracy')

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

if __name__ == "__main__":
	(options, args) = argsParse()
	In = options.input
	script = options.script
	Out_tmp = options.out + ".tmp"
	length = options.length
	degeneracy = options.deg
	e1 = time.time()
	os.system("perl {}/DEGEPRIME-1.1.0/DegePrime.pl -i {} \
                        -d {} -l {} -o {}".format(script,In,degeneracy,length,Out_tmp))
	
	os.rename(Out_tmp, options.out)
	e2 = time.time()
	print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))















