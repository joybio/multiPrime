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

parser.add_option('-i','--input',
                dest='input',
                help='input file: primer.fa.')

parser.add_option('-r','--ref',
                dest='ref',
                help='reference file: bowtie2.')

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
import os
import threading
import time
import pandas as pd
from time import strftime
import functools
def get_term9(fasta,out):
	with open(fasta,"r") as f:
		with open(out,"w") as o:
			for i in f:
				if i.startswith(">"):
					o.write(i)
				else:
					i = i.strip()[-9:]
					o.write(i + "\n")
def map(fasta,path_2_ref,out,for_out,rev_out):
	os.system("bowtie2 -p 20 -f -N 0 -a -x {} -f -U {} -S {}".format(path_2_ref,fasta, out))
	os.system("samtools view -F 16 {} > {}".format(out,for_out))
	os.system("samtools view -f 16 {} > {}".format(out,rev_out))

def build_dict(input_file, input_dict,strand):
	global prediction	
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
								prediction = pd.concat([f_prediction,offtarget],axis=0,ignore_index=True)
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
								prediction = pd.concat([r_prediction,offtarget],axis=0,ignore_index=True)
if __name__ == '__main__':
	print("INFO {} Start: Load file ...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
	product_len = options.product.split(",")
	print("offtargets length: {} - {}".format(product_len[0],product_len[1]))
	fasta = options.input
	term9_fasta = options.input.rstrip("fa")+ "9.fa"
	get_term9(fasta,term9_fasta)
	path_2_ref = options.ref
	sam_file = options.input.rstrip("fa")+ "sam"
	sam_for_file = options.input.rstrip("fa")+ "for.sam"
	sam_rev_file = options.input.rstrip("fa")+ "rev.sam"
	map(term9_fasta,path_2_ref,sam_file,sam_for_file,sam_rev_file)
	forward_dict = defaultdict(list)
	reverse_dict = defaultdict(list)
	prediction = pd.DataFrame(columns=["Primer_F","Primer_R","Product length","Chrom (or Genes)","Start","Stop"])
	t = []
	t.append(threading.Thread(target=build_dict, args=(sam_for_file, forward_dict,"F")))
	t.append(threading.Thread(target=build_dict, args=(sam_rev_file, reverse_dict,"R")))
	for t1 in t:
		t1.start()
	for t1 in t:
		t1.join()
	with open(options.out,"w") as f:
		prediction.to_csv(f, index=False, sep="\t")
	print("INFO {} Done...".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))


