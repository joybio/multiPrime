#!/bin/pyhthon
"""genome format"""
#caution the chrom column of the position file must contain chr
__date__ = "2022-6-1"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input.fa] -o [output.fa]')
parser.add_option('-i','--input',
		dest='input',
		help='Input file: Check')

parser.add_option('-f','--format',
		dest='format',
		help='Format of primer file: xls or fa or seq; default: xls. \n xls: final_primer_set.xls, output of multiPrime. \n \
			 fa: fasta format. \n seq: sequence format, comma seperate. e.g. primer_F,Primer_R.')

parser.add_option('-p','--primer',
		dest='primer',
		help='Primer file. One of the followed three types:\n final_maxprimers_set.xls \n primer.fa \n primer_F,primer_R.')

parser.add_option('-o','--out',
		dest='out',
		help='Out file: fasta')

parser.add_option('-s','--stast',
		dest='stast',
		default="Coverage.xls",
		help='Stast information: number of coverage and total. default: Coverage.xls')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
import os 
os.system("cat {} | grep {} > {}".format(options.input,,))
os.system("bedtools getfasta -fi {} -bed {} -fo {}".format())
def get_PCR_PRODUCT(ref,F,R):
	Fseq = dege_trans(F)
	Rseq = dege_trans(R)
	product_dict ={}
	global target



if __name__ == "__main__":
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


