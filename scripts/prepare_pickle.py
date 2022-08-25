#!/bin/python
"""
build a binary dict using the columns of input file
"""

__date__ = "2022-6-1"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import pickle
import optparse
import collections
from collections import defaultdict
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input] -n [int] -o [output] options {-v [int or T]} \n return a binary dict, which can loaded by pickle')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: List of virus')

parser.add_option('-n','--num',
                dest='num',
                help='number: which column should be used as hash key. for example: \
			if you want to set column 1 as the key of dict, then -n 0.')

parser.add_option('-v','--value',
                dest='value',
		default="T",
                help='value: which column should be used as hash value.\
			if value == T(total),use the whole line as the hash value')

parser.add_option('-o','--out',
                dest='out',
                help='Out file: virus genome.fa')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)


data = open(options.input,"r")
out = open(options.out,"wb")
column = int(options.num)
value = options.value
print(value)
pickle_dict = defaultdict(list)
for i in data:
	line = i.strip()
	i = i.strip().split("\t")
	#print(i)
	key = i[column]
	if value == "T":
		#print(key,line)
		pickle_dict[key]=line
	else:
		pickle_dict[key].append(i[int(value)])
		#print(key,i[int(value)])
pickle.dump(pickle_dict,out)
data.close()
out.close()

