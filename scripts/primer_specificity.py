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

parser = OptionParser('Usage: %prog -i [input] -o [output]')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: blast out.')
parser.add_option('-t','--term',
                dest='term',
		default="0",
		type="int",
                help='distance to term, default: 0')
parser.add_option('-l','--length',
                dest='length',
                default="8",
                type="int",
                help='match length, default: 8')
parser.add_option('-m','--mis',
                dest='mismatch',
                default="1",
                type="int",
                help='mismatch number, default: 1')
parser.add_option('-p','--product',
                dest='product',
                default="20,2000",
                help='length of product, default: [20,2000]')
parser.add_option('-o','--out',
                dest='out',
                help='non-specifict PCR product')

(options,args) = parser.parse_args()
if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
out = open(options.out,"w")
import collections
from collections import defaultdict
import re
dist_term = options.term
length = options.length
gene_start_dict = defaultdict(list)
gene_stop_dict = defaultdict(list)
mismatch = options.mismatch
product_len = options.product.split(",")
with open(options.input,"r") as f:
	for i in f:
		i = i.strip().split("\t")
		gene = i[4]
		primer = i[0]
		mismatch = int(i[14])
		primer_start = i[2]
		primer_stop = i[3]
		primer_len = i[1]
		if int(i[1])-int(i[3]) <= dist_term and int(i[3]) - int(i[2]) > length and mismatch <= options.mismatch:
			if i[6] < i[7]:
				match_start = i[6]
				match_stop = i[7]
				gene_start_dict[gene].append([primer,match_start,primer_start,primer_stop,primer_len,match_stop,mismatch])
			else:
				match_start = i[7]
				match_stop = i[6]
				gene_stop_dict[gene].append([primer,match_start,primer_start,primer_stop,primer_len,match_stop,mismatch])
out.write("Chrom\tProduct_start\tProduct_stop\tProduct_length\tPrimer_F\tPrimer_F_len\tPrimer_F_start\tPrimer_F_stop\tPrimer_F_mis\tPrimer_R\tPrimer_R_len\tPrimer_R_start\tPrimer_R_stop\tPrimer_R_mis\tPrimerF_start_pos\tPrimerF_stop_pos\tPrimerR_start_pos\tPrimerR_stop_pos\n")
for i in gene_start_dict.keys():
	if i in gene_stop_dict.keys():
		for start in gene_start_dict[i]:
			primer_F = start[0]
			primerF_start_pos = int(start[1])
			primerF_mis_start_stop = start[2] + "\t" + start[3]
			primer_F_len = start[4]
			primerF_stop_pos = start[5]
			primerF_mismatch = start[6]
			for stop in gene_stop_dict[i]:
				primer_R = stop[0]
				primerR_stop_pos = int(stop[1])
				primerR_mis_start_stop = stop[2] + "\t" + stop[3]
				primer_R_len = stop[4]
				primerR_start_pos = int(stop[5])
				primerR_mismatch = start[6]
				if primerR_start_pos - primerF_start_pos > int(product_len[0]) and primerR_start_pos - primerF_start_pos < int(product_len[1]):
					product_start = str(primerF_start_pos)
					product_stop = str(primerR_start_pos)
					out.write(i + "\t" + product_start + "\t" + product_stop +"\t" + str(primerR_start_pos-primerF_start_pos) +"\t" + \
						primer_F +"\t" + primer_F_len +"\t" + primerF_mis_start_stop+"\t"+ str(primerF_mismatch) +"\t"+ \
						primer_R + "\t" + primer_R_len +"\t" + primerR_mis_start_stop + "\t" + str(primerR_mismatch) + "\t" + \
						str(primerF_start_pos) + "\t" + str(primerF_stop_pos) + "\t" + str(primerR_start_pos) + "\t" + str(primerR_stop_pos) + "\n")
out.close()		
