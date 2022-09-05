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
import pandas as pd
from time import strftime
import functools
product_len = options.product.split(",")

def build_dict(input_file, input_dict,strand):
	global f_prediction, r_prediction	
	with open(input_file,"r") as f:
		for i in f:
			i = i.strip().split("\t")
			primer = i[0]
			gene = i[2]
			primer_match_start = int(i[3])-1
			input_dict[gene].append([primer,primer_match_start])
	if strand == "F":
		for gene in input_dict.keys():
			for primer_F in input_dict[gene]:
				if primer_F[0].endswith("_R"):
					#print(primer_F)
					position = input_dict[gene].index(primer_F)
					for primer_R in input_dict[gene][position+1:]:
						if primer_R[0].endswith("_F"):
							distance = int(primer_R[1])-int(primer_F[1])+1
							if distance > int(product_len[0]) and distance < int(product_len[1]):
								offtarget = pd.DataFrame({"Primer_F":primer_F[0],
										"Primer_R":primer_R[0],
										"Product length":distance,
										"Chrom (or Genes)":gene,
										"Start":int(primer_F[1])-1,
										"Stop":int(primer_R[1])-1
										},index=[0])
								f_prediction = pd.concat([f_prediction,offtarget],axis=0,ignore_index=True)
	elif strand == "R":
		for gene in input_dict.keys():
			for primer_F in input_dict[gene]:
				if primer_F[0].endswith("_F"):
					position = input_dict[gene].index(primer_F)
					for primer_R in input_dict[gene][position+1:]:
						if primer_R[0].endswith("_R"):
							distance = int(primer_R[1])-int(primer_F[1])+1
							if distance > int(product_len[0]) and distance < int(product_len[1]):
								offtarget = pd.DataFrame({"Primer_F":primer_F[0],
									"Primer_R":primer_R[0],
									"Product length":distance,
									"Chrom (or Genes)":gene,
									"Start":int(primer_F[1])-1,
									"Stop":int(primer_R[1])-1
									},index=[0])
								r_prediction = pd.concat([r_prediction,offtarget],axis=0,ignore_index=True)
if __name__ == '__main__':
	print("INFO {} Start: Load file ...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
	forward_dict = defaultdict(list)
	reverse_dict = defaultdict(list)
	f_prediction = pd.DataFrame(columns=["Primer_F","Primer_R","Product length","Chrom (or Genes)","Start","Stop"])
	r_prediction = pd.DataFrame(columns=["Primer_F","Primer_R","Product length","Chrom (or Genes)","Start","Stop"])
	#global f_prediction, r_prediction	
	t = []
	t.append(threading.Thread(target=build_dict, args=(options.forward, forward_dict,"F")))
	t.append(threading.Thread(target=build_dict, args=(options.reverse, reverse_dict,"R")))
	for t1 in t:
		t1.start()
	for t1 in t:
		t1.join()
	print(f_prediction)
	prediction = pd.concat([f_prediction,r_prediction], axis=0, ignore_index=True)
	with open(options.out,"w") as f:
		prediction.to_csv(f, index=False, sep="\t")
	print("INFO {} Done...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))


