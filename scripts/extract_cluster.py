#!/bin/python

"""extract clusters from ch-hit results"""
# caution the chrom column of the position file must contain chr
__date__ = "2022-6-7"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
from optparse import OptionParser
import random
from collections import defaultdict
import re
import os


def argsParse():
	parser = OptionParser('Usage: %prog -i [total.fa] -c [cluster.clstr] \n'
						  'Options: -o [cluster.txt] -y [cluster.identities] -m [500] -d [Cluster.fa]')
	parser.add_option('-i', '--input',
					  dest='input',
					  help='Input file: *complete.genome.fa')
	parser.add_option('-c', '--cluster',
					  dest='cluster',
					  help='cluster file: *.clstr')
	parser.add_option('-o', '--out',
					  dest='out',
					  default="cluster.txt",
					  type="str",
					  help='output file: clusters information. Default: cluster.txt. \
						It will be used as input file of snakemake pipeline')
	parser.add_option('-y', '--identity',
					  dest='identity',
					  default="cluster.identities.txt",
					  type="str",
					  help='output file: clusters information. Default: cluster.identities.txt. \
						It will be used as input file of snakemake pipeline')
	parser.add_option('-m', '--max',
					  dest='max',
					  default=500,
					  type="int",
					  help='max sequence number in 1 cluster. Default: 500.')
	parser.add_option('-d', '--dir',
					  dest='dir',
					  default="Cluster_fa",
					  type="str",
					  help='directory of output fasta: clusters information. Default: Cluster_fa.')
	(options, args) = parser.parse_args()
	import sys
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	elif options.input is None:
		parser.print_help()
		print("Input file must be specified !!!")
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

if __name__ == "__main__":
	(options, args) = argsParse()
	data = open(options.input, "r")
	cluster = open(options.cluster, "r")
	out = open(options.out, "w")
	directory = options.dir
	cluster_dict = defaultdict(list)
	cluster_identity = defaultdict(list)
	cluster_rep = {}
	seq_dict = {}
	identities = open(options.identity, "w")
	for c in cluster:
		if c.startswith(">"):
			cluster_name = c.strip().replace(" ", "_")
		else:
			c = c.strip().split(">")
			value = c[1].split("...")
			cluster_dict[cluster_name].append(">" + value[0])
			identity = c[1].split(" ")
			cluster_identity[cluster_name].append(identity[-1])
			if identity[-1] == "*":
				cluster_rep[cluster_name] = ">" + value[0]
	cluster.close()
	for i in cluster_identity.keys():
		for j in cluster_identity[i]:
			if j != "*":
				identities.write(i.lstrip(">") + '\t' + j + "\n")
	identities.close()

	for seq in data:
		if seq.startswith(">"):
			if re.search(" ", seq):
				seq = seq.split(" ")
				seq_name = seq[0]
			else:
				seq_name = seq.strip()
		else:
			seq_dict[seq_name] = seq
	data.close()
	for k in cluster_dict.keys():
		os.system("mkdir -p {}".format(directory))
		out.write(k.lstrip(">") + "_" + str(len(cluster_dict[k])) + ".fa\n")
		# out.write(directory + "/" + k.lstrip(">")+"_"+str(len(cluster_dict[k])) + "/" +k.lstrip(">")+"_"+str(len(cluster_dict[k]))+".fa\n")
		cluster_id = directory + "/" + k.lstrip(">") + "_" + str(len(cluster_dict[k])) + ".fa"
		with open(cluster_id, "w") as f:
			for l in cluster_dict[k]:
				f.write(str(l) + "\n" + seq_dict[l])
		top_cluster_id = directory + "/" + k.lstrip(">") + "_" + str(len(cluster_dict[k])) + ".tfa"
		if len(cluster_dict[k]) < int(options.max):
			with open(top_cluster_id, "w") as t:
				for l in cluster_dict[k]:
					t.write(str(l) + "\n" + seq_dict[l])
		else:
			cluster_dict[k].remove(cluster_rep[k])
			selected_seq = random.sample(cluster_dict[k], int(options.max) - 1)
			selected_seq.append(cluster_rep[k])
			with open(top_cluster_id, "w") as t:
				for seq_id in set(selected_seq):
					t.write(str(seq_id) + "\n" + seq_dict[seq_id])
		cluster_seq = directory + "/" + k.lstrip(">") + "_" + str(len(cluster_dict[k])) + ".txt"
		with open(cluster_seq, "w") as s:
			for i in cluster_dict[k]:
				s.write(i + "\n")

	out.close()
