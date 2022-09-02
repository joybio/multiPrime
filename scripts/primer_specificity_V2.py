#!/bin/python
"""
Output off-target PCR results. 
"""
__date__ = "2022-8-15"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import sys
from sys import argv
import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input] -f [forward.sam] -r [reverse.sam] -p [50,2000]-o [output]')

parser.add_option('-f','--for',
                dest='forward',
                help='forward file: forward.sam.')

parser.add_option('-r','--rev',
                dest='reverse',
                help='reverse file: reverse.sam.')

parser.add_option('-p','--product',
                dest='product',
                default="50,2000",
                help='length of product, default: [50,2000]')
parser.add_option('-o','--out',
                dest='out',
                help='non-specifict PCR product')

(options,args) = parser.parse_args()
if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
import collections
from collections import defaultdict
import re
import threading
import time
from time import strftime
import functools
product_len = options.product.split(",")

def build_dict(input_file, input_dict):
	with open(input_file,"r") as f:
		for i in f:
			i = i.strip().split("\t")
			primer = i[0]
			gene = i[2]
			primer_match_start = int(i[3])-1
			input_dict[gene].append([primer,primer_match_start])

if __name__ == '__main__':
	print("INFO {} Start: Load file ...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
	out = open(options.out,"w")
	forward_dict = defaultdict(list)
	reverse_dict = defaultdict(list)
	print(options.forward)
	t = []
	t.append(threading.Thread(target=build_dict, args=(options.forward, forward_dict)))
	t.append(threading.Thread(target=build_dict, args=(options.reverse, reverse_dict)))
	for t1 in t:
		t1.start()
	for t1 in t:
		t1.join()
	print(len(forward_dict))
	out.write("offtarget_primers\tproduct_len\tofftarget_chrom\tofftarget_start\tofftarget_stop\n")
	print("INFO {} Loding finished ...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
	for target in forward_dict.keys():
		if target in reverse_dict:
			for primer_start in forward_dict[target]:
				for primer_stop in reverse_dict[target]:
					target_start = primer_start[1]
					target_stop = primer_stop[1]
					if target_stop - target_start < 2000 and target_stop - target_start > 50:
						primer_F = primer_start[0]
						primer_R = primer_stop[0]
						product_len = target_stop - target_start
						out.write(primer_F+":"+primer_R+"\t"+str(product_len) + "\t"+target+"\t"+str(target_start)+"-"+str(target_stop)+"\n")
	out.close()
	print("INFO {} Done...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))


