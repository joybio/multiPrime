#!/bin/python
'''extract PCR product by primers (degenerate is also ok!) and raw input fasta file.'''
import re
import sys
from sys import argv
import optparse
from optparse import OptionParser
from Bio.Seq import Seq
parser = OptionParser('Usage: %prog -i [input] -p [primerF,primerR] -o [output]')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: ref_fa.')

parser.add_option('-p','--primer',
                dest='primer',
                help='primers: final_maxprimers_set.xls')

parser.add_option('-o','--out',
                dest='out',
                help='output_dir')

parser.add_option('-s','--stast',
                dest='stast',
                help='stast information: coverage/total')

(options,args) = parser.parse_args()

# replace degenerate base to [A,C,G,T].
import math
import os
from operator import mul
from functools import reduce

degenerate_pair = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

degenerate_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
                    "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}

'''
def get_dict_key(d_dict, value):
    d_key = list(d_dict.keys())[list(d_dict.values()).index(value)]  # 如果两个key对应同一个value，那必然报错。
    return d_key
'''


def score_trans(sequence):
    return reduce(mul, [math.floor(degenerate_table[x]) for x in list(sequence)])


def dege_trans(sequence):
    expand_seq = [sequence]
    expand_score = reduce(mul, [score_trans(x) for x in expand_seq])
    while expand_score > 1:
        for seq in expand_seq:
            if score_trans(seq) == 1:
                pass
            else:
                expand_seq.remove(seq)
                for nt in range(len(seq)):
                    if math.floor(degenerate_table[seq[nt]]) > 1:
                        for v in degenerate_pair[seq[nt]]:
                            expand_seq.append(seq[0:nt] + v + seq[nt + 1:])
                            # print(expand_seq)
                        break
        expand_score = reduce(mul, [score_trans(x) for x in expand_seq])
    return expand_seq

if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
primer = open(options.primer,"r")
reference = options.input

def get_PCR_PRODUCT(ref,F,R):
	Fseq = dege_trans(F)
	Rseq = dege_trans(R)
	product_dict ={}
	with open(ref,"r") as r:
		for i in r:
			if i.startswith(">"):
				key = i.strip()
			else:
				for sequence in Fseq:
					value = ''
					if re.search(sequence,i):
						line = i.split(sequence)
						product = sequence+line[1]
						for sequence2 in Rseq:
							if re.search(str(Seq(sequence2).reverse_complement()),product):
								product = product.split(str(Seq(sequence2).reverse_complement()))
								value = product[0].strip() + str(Seq(sequence2).reverse_complement())
								product_dict[key] = value
								break
							else:
								pass
						if value:
							break
					else:
						pass

	return product_dict

if __name__ == "__main__":
	seq_id = set()
	coverage_number = open(options.stast,"w")
	Total_seq_number = 0
	for i in primer:
		if i.startswith("#"):
			pass
		else:
			i = i.strip().split("\t")
			cluster_id = i[0].split("/")
			cluster_number = cluster_id[-1].rstrip(".candidate.primers.txt").split("_")
			Total_seq_number += int(cluster_number[-1])
			out_dir = options.out
			os.system("mkdir -p {}".format(out_dir))
			cluster_product = options.out + "/" + cluster_id[-1].rstrip(".candidate.primers.txt") + ".PCR.priduct.fa"
			primer_F = i[2]
			primer_R = i[3]
			PCR_product = get_PCR_PRODUCT(reference,primer_F,primer_R)
			with open(cluster_product,"w") as c:
				for i in PCR_product.keys():
					c.write(i + "\n" + PCR_product[i] + "\n")
					seq_id.add(i)
	coverage_number.write("Total number of sequences:\t{}\nCoveraged number of sequence\t{}\nRate of coverage:\t{}\n".format(Total_seq_number,len(seq_id),round(float(len(seq_id))/Total_seq_number,2)))
	coverage_number.close()	






