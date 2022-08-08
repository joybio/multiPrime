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
                help='Input file: degeprimer out.')

parser.add_option('-p','--primer',
                dest='primer',
                help='primer sequence: seperate by comma')
parser.add_option('-o','--out',
                dest='out',
                help='Prefix of out file: candidate primers')

(options,args) = parser.parse_args()

# replace degenerate base to [A,C,G,T].
import math
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
raw_fa = open(options.input,"r")
primer = options.primer.split(",")
Fseq = dege_trans(primer[0])
Rseq = dege_trans(primer[1])

product_dict ={}
for i in raw_fa:
	if i.startswith(">"):
		key = i.strip()
	else:
		for sequence in Fseq:
			value = ''
			#print(sequence,i)
			if re.search(sequence,i):
				line = i.split(sequence)
				product = sequence+line[1]
				for sequence2 in Rseq:
					#print(str(Seq(sequence2).reverse_complement()),product)
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
out = open(options.out,"w")
for i in product_dict.keys():
	out.write(i + "\n" + product_dict[i] + "\n")	
out.close()






