#!/share/data3/yangjunbo/miniconda3/bin/python
'''extract perfect PCR product by primers (degenerate is also ok!) with raw input fasta file.'''
import re
import sys
from sys import argv
import optparse
from optparse import OptionParser
from Bio.Seq import Seq
parser = OptionParser('Usage: %prog -i [input] -p [primerF,primerR] -f [format] -o [output]',version = "%prog 0.0.2")
parser.add_option('-i','--input',
                dest='input',
                help='Input file: template fasta or reference fasta.')

parser.add_option('-p','--primer',
                dest='primer',
                help='Primer file. One of the followed three types:\n final_maxprimers_set.xls \n primer.fa \n primer_F,primer_R.')

parser.add_option('-f','--format',
                dest='format',
                help='Format of primer file: xls or fa or seq; default: xls. \n xls: final_primer_set.xls, output of multiPrime. \n \
			 fa: fasta format. \n seq: sequence format, comma seperate. e.g. primer_F,Primer_R.')

parser.add_option('-o','--out',
                dest='out',
		default="PCR_product",
                help='Output_dir. default: PCR_product.')

parser.add_option('-s','--stast',
                dest='stast',
                default="Coverage.xls",
		help='Stast information: number of coverage and total. default: Coverage.xls')
(options,args) = parser.parse_args()

# replace degenerate base to [A,C,G,T].
import math
import os
from operator import mul
from functools import reduce
import pandas as pd

degenerate_pair = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

degenerate_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
                    "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}

'''
def get_dict_key(d_dict, value):
    d_key = list(d_dict.keys())[list(d_dict.values()).index(value)]  
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

def get_PCR_PRODUCT(ref,F,R):
	Fseq = dege_trans(F)
	Rseq = dege_trans(R)
	product_dict ={}
	global target
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
								target += 1
								break
							else:
								pass
						if value:
							break
					else:
						pass

	return product_dict

if __name__ == "__main__":
	target,total = 0,0
	seq_id = set()
	reference = options.input
	out_dir = options.out
	os.system("mkdir -p {}".format(out_dir))
	coverage_number = open(options.stast,"w")
	if options.format == 'seq':
		seq_product = options.out + "/" + "PCR.product.fa"
		primers = options.primer.split(",")
		primer_F = primers[0]
		primer_R = primers[1]
		PCR_product = get_PCR_PRODUCT(reference,primer_F,primer_R)
		with open(seq_product,"w") as s:
			for i in PCR_product.keys():
				s.write(i + "\n" + PCR_product[i] + "\n")
				seq_id.add(i)
			coverage_number.write("Number of {}:\t{}\n".format(options.primer,target))
			target = 0
	else:
		primer = open(options.primer,"r")
		if options.format == 'xls':
			coverage_number = open(options.stast,"w")
			for i in primer:
				if i.startswith("#"):
					pass
				else:
					i = i.strip().split("\t")
					cluster_id = i[0].split("/")		
					cluster_product = options.out + "/" + cluster_id[-1].rstrip(".candidate.primers.txt") + ".PCR.product.fa"
					primer_F = i[2]
					primer_R = i[3]
					PCR_product = get_PCR_PRODUCT(reference,primer_F,primer_R)
					with open(cluster_product,"w") as c:
						for p in PCR_product.keys():
							c.write(p + "\n" + PCR_product[p] + "\n")
							seq_id.add(p)
					coverage_number.write("Number of {}:\t{}\n".format(cluster_product,target))
					target = 0
		elif options.format == 'fa':
			primer_info = pd.read_table(primer,header=None)
			for idx,row in primer_info.iterrows():
				if idx%4==0:
					primer_F_info = options.out + "/" + row[0].lstrip(">")
				elif idx%4==1:
					primer_F = row[0]
				elif idx%4==2:
					fa_product = primer_F_info + "_" + row[0].lstrip(">")+ ".PCR.product.fa"
				elif idx%4==3:
					primer_R = row[0]
					PCR_product = get_PCR_PRODUCT(reference,primer_F,primer_R)
					with open(fa_product,"w") as f:
						for i in PCR_product.keys():
							f.write(i + "\n" + PCR_product[i] + "\n")
							seq_id.add(i)
					coverage_number.write("Number of {}:\t{}\n".format(fa_product,target))
					target = 0
		else:
			parser.print_help()
			print("Please check you primer format!")
			sys.exit(1)
	with open(reference,"r") as r:
		for i in r:
			if i.startswith(">"):
				total += 1
	coverage_number.write("Total number of sequences:\t{}\nCoveraged number of sequence:\t{}\nRate of coverage:\t{}\n".format(total,len(seq_id),round(float(len(seq_id))/total,2)))
	coverage_number.close()



