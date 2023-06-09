#!/share/data3/yangjunbo/miniconda3/bin/python

"""
extract value from a subset of key
"""
__date__ = "2022-6-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
import pickle
from optparse import OptionParser
import os
import re
import time
from time import strftime

parser = OptionParser('Usage: %prog  -i [input] -v [virus_dict] -o [output] \n \
		Options: -n [num] \n \
		extract value from a subset of key. It accept duplicated keys[means more than 1 key in column], values will be a list,separated by line break!')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: a subset of key, need values.')

parser.add_option('-n','--num',
                dest='num',
		default=0,
                help='Column number: which Colum is the keys. Default: [0],first column')

parser.add_option('-d','--dict',
                dest='dict',
                default = "/share/data3/yangjunbo/database/accession_taxid.dict",
                help='defualt dict: /share/data3/yangjunbo/database/accession_taxid.dict; virus hash: keys = acc; values = line. \n \
			you can construct your own dict by prepare_pickle.py.')

parser.add_option('-o','--out',
                dest='out',
                help='Out file')
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

(options,args) = parser.parse_args()
dictionary = open(options.dict,"rb")
raw_dict = pickle.load(dictionary)
number = int(options.num)
data = open(options.input,"r")
out = open(options.out,"w")
print(raw_dict)
for i in data:
	i = i.strip().split("\t")
	key = i[number]
	if key in raw_dict.keys():
		#print(raw_dict[key])
		out.write(raw_dict[key] + "\n")
	else:
		pass
data.close()
dictionary.close()
out.close()




