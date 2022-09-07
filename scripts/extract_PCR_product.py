#!/share/data3/yangjunbo/miniconda3/bin/python
'''extract perfect PCR product by primers (degenerate is also ok!) with raw input fasta file.'''
import re
import sys
from optparse import OptionParser
from Bio.Seq import Seq
from get_degePrimer import dege_trans
import os
from operator import mul
from functools import reduce
import pandas as pd


def argsParse():
	parser = OptionParser('Usage: %prog -i [input] -p [primerF,primerR] -f [format] -o [output]', version="%prog 0.0.2")
	parser.add_option('-i', '--input',
					  dest='input',
					  help='Input file: template fasta or reference fasta.')

	parser.add_option('-p', '--primer',
					  dest='primer',
					  help='Primer file. One of the followed three types:\n final_maxprimers_set.xls \n primer.fa \n primer_F,primer_R.')

	parser.add_option('-f', '--format',
					  dest='format',
					  help='Format of primer file: xls or fa or seq; default: xls. \n xls: final_primer_set.xls, output of multiPrime. \n \
				 fa: fasta format. \n seq: sequence format, comma seperate. e.g. primer_F,Primer_R.')

	parser.add_option('-o', '--out',
					  dest='out',
					  default="PCR_product",
					  help='Output_dir. default: PCR_product.')

	parser.add_option('-s', '--stast',
					  dest='stast',
					  default="Coverage.xls",
					  help='Stast information: number of coverage and total. default: Coverage.xls')
	(options, args) = parser.parse_args()
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	elif options.input is None:
		parser.print_help()
		print("Input file must be specified !!!")
		sys.exit(1)
	elif options.primer is None:
		parser.print_help()
		print("Primer file or sequence must be specified !!!")
		sys.exit(1)
	elif options.format is None:
		parser.print_help()
		print("Primer file format must be specified !!!")
		sys.exit(1)
	elif options.out is None:
		parser.print_help()
		print("No output file provided !!!")
		sys.exit(1)
	count = 3  # times for the packages install
	while count:
		try:
			import Bio  #

			print('Dependent package Biopython is OK.\nDependent module Bio is OK.')
			break
		except:
			print('Dependent package Biopython is not found!!! \n Start installing ....')
			os.system('pip install biopython')
			count -= 1
			continue
	return parser.parse_args()


def get_PCR_PRODUCT(ref, F, R):
	Fseq = dege_trans(F)
	Rseq = dege_trans(R)
	product_dict = {}
	global target
	with open(ref, "r") as r:
		for i in r:
			if i.startswith(">"):
				key = i.strip()
			else:
				for sequence in Fseq:
					value = ''
					if re.search(sequence, i):
						line = i.split(sequence)
						product = sequence + line[1]
						for sequence2 in Rseq:
							if re.search(str(Seq(sequence2).reverse_complement()), product):
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
	(options, args) = argsParse()
	target, total = 0, 0
	seq_id = set()
	reference = options.input
	out_dir = options.out
	os.system("mkdir -p {}".format(out_dir))
	coverage_number = open(options.stast, "w")
	if options.format == 'seq':
		seq_product = options.out + "/" + "PCR.product.fa"
		primers = options.primer.split(",")
		primer_F = primers[0]
		primer_R = primers[1]
		PCR_product = get_PCR_PRODUCT(reference, primer_F, primer_R)
		with open(seq_product, "w") as s:
			for i in PCR_product.keys():
				s.write(i + "\n" + PCR_product[i] + "\n")
				seq_id.add(i)
			coverage_number.write("Number of {}:\t{}\n".format(options.primer, target))
			target = 0
	else:
		primer = open(options.primer, "r")
		if options.format == 'xls':
			coverage_number = open(options.stast, "w")
			for i in primer:
				if i.startswith("#"):
					pass
				else:
					i = i.strip().split("\t")
					cluster_id = i[0].split("/")
					cluster_product = options.out + "/" + cluster_id[-1].rstrip(
						".candidate.primers.txt") + ".PCR.product.fa"
					primer_F = i[2]
					primer_R = i[3]
					PCR_product = get_PCR_PRODUCT(reference, primer_F, primer_R)
					with open(cluster_product, "w") as c:
						for p in PCR_product.keys():
							c.write(p + "\n" + PCR_product[p] + "\n")
							seq_id.add(p)
					coverage_number.write("Number of {}:\t{}\n".format(cluster_product, target))
					target = 0
		elif options.format == 'fa':
			primer_info = pd.read_table(primer, header=None)
			for idx, row in primer_info.iterrows():
				if idx % 4 == 0:
					primer_F_info = options.out + "/" + row[0].lstrip(">")
				elif idx % 4 == 1:
					primer_F = row[0]
				elif idx % 4 == 2:
					fa_product = primer_F_info + "_" + row[0].lstrip(">") + ".PCR.product.fa"
				elif idx % 4 == 3:
					primer_R = row[0]
					PCR_product = get_PCR_PRODUCT(reference, primer_F, primer_R)
					with open(fa_product, "w") as f:
						for i in PCR_product.keys():
							f.write(i + "\n" + PCR_product[i] + "\n")
							seq_id.add(i)
					coverage_number.write("Number of {}:\t{}\n".format(fa_product, target))
					target = 0
		else:
			print("Please check you primer format!")
			sys.exit(1)
	with open(reference, "r") as r:
		for i in r:
			if i.startswith(">"):
				total += 1
	coverage_number.write(
		"Total number of sequences:\t{}\nCoveraged number of sequence:\t{}\nRate of coverage:\t{}\n".format(total,
																											len(seq_id),
																											round(float(
																												len(seq_id)) / total,
																												  2)))
	coverage_number.close()
