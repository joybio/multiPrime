#!/bin/python
"""
Selcet non-specific PCR results.
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
		default="1",
		type="int",
                help='distance to term, default: 1')
parser.add_option('-l','--length',
                dest='length',
                default="8",
                type="int",
                help='match length, default: 8')
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
with open(options.input,"r") as f:
	for i in f:
		i = i.strip().split("\t")
		gene = i[4]
		primer = i[0]
		mismatch = int(i[14])
		if int(i[1])-int(i[3]) <= dist_term and int(i[3]) - int(i[2]) > length and mismatch <= 1:
			if i[6] < i[7]:
				primer_start = i[6]
				gene_start_dict[gene].append([primer,primer_start,i[2],i[3],i[1],i[7]])
			else:
				primer_start = i[7]
				gene_stop_dict[gene].append([primer,primer_start,i[2],i[3],i[1],i[6]])

out.write("Chrom\tPrimer_F_start\tPrimer_R_stop\tProduct_length\tPrimer_F\tPrimer_F_len\tPrimer_F_start\tPrimer_F_stop\tPrimer_R\tPrimer_R_len\tPrimer_R_start\tPrimer_R_stop\n")
for i in gene_start_dict.keys():
	if i in gene_stop_dict.keys():
		for start in gene_start_dict[i]:
			primer_F = start[0]
			start_pos = int(start[1])
			primerF_mis_start_stop = start[2] + "\t" + start[3]
			primer_F_len = start[4]
			start_pos2 = start[5]
			for stop in gene_stop_dict[i]:
				primer_R = stop[0]
				stop_pos = int(stop[1])
				primerR_mis_start_stop = start[2] + "\t" + start[3]
				primer_R_len = stop[4]
				stop_pos2 = int(stop[5])
				if re.search("_R",primer_F) and re.search("_R",primer_R):
					pass
				else:
					if stop_pos - start_pos > -2000 and stop_pos -start_pos < -20:
						product_start = str(stop_pos)
						product_stop = str(start_pos2)
						out.write(i + "\t" + product_start + "\t" + product_stop +"\t" + str(abs(stop_pos2-start_pos)) +"\t" + \
							primer_F +"\t" + primer_F_len +"\t" + primerF_mis_start_stop +"\t"+ \
							primer_R + "\t" + primer_R_len +"\t" + primerR_mis_start_stop + "\t" + str(start_pos) + \
							"\t" + str(start_pos2) + "\t" + str(stop_pos2) + "\t" +  str(stop_pos) + "\n")
					elif stop_pos - start_pos > 20 and stop_pos -start_pos < 2000:
						product_start = str(start_pos)
						product_stop = str(stop_pos2)
						out.write(i + "\t" + product_start + "\t" + product_stop +"\t" + str(abs(stop_pos2-start_pos)) +"\t" + \
							primer_F +"\t" + primer_F_len +"\t" + primerF_mis_start_stop +"\t"+ \
							primer_R + "\t" + primer_R_len +"\t" + primerR_mis_start_stop + "\t" + str(start_pos) + \
							"\t" + str(start_pos2) + "\t" + str(stop_pos2) + "\t" +  str(stop_pos) + "\n")
out.close()		
