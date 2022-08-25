#!/bin/python
"""
Selcet candidate primers from degePrimer results.
"""
__date__ = "2022-6-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import sys
from sys import argv
import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input] -r [sequence.fa] -o [output] \n \
			Options: {-f [fraction] -n [num; -f 0] -s [150,500] -g [0.45,0.65] -d [ditance] -k [kmer] -w [window]}. \n \
			If there are too many sequence in the reference fasta. You can use the first [N;N<=50,corresponding to the param -n] sequence for the degeprimer input. \n  \
			use degePrimer output as -i [degeprimer_output] and all sequence as -r [all.sequence.fa]. get_degePrimer.py will help you to choose the candidate primers.')
parser.add_option('-i','--input',
		dest='input',
		help='Input file: degeprimer out.')

parser.add_option('-r','--ref',
		dest='ref',
		help='Reference sequence file: all the sequence in 1 fasta, for example: (Cluster_96_171.fa)')

parser.add_option('-f','--frac',
		dest='frac',
		default = "0.9",
		help='Filter primers by match rate: [Number of sequences that match the selected primer] / [Number of sequences that span the selected primer]. Filter primer by fraction. Only primers with fraction [> frac] will retain. default: 0.9.')

parser.add_option('-n','--num',
		dest='num',
		default = "1",
		help='Filter primers by match number: Number of sequences that match the selected primer. \n \
		Sometimes, the input sequence is too many, and a subset of the total input sequences is used. Under this condition, if you use the fraction param, It wont return candidate primers, because the match rate is very small, even nearly zero. However, you can use the number of subset sequence to select candidate primers. \n \
		before use this param, please also use [ -f 0 ] \n \
		default: 1')

parser.add_option('-g','--gc',
		dest='gc',
		default = "0.45,0.65",
		help="Filter primers by GC content. default [0.45,0.65].")

parser.add_option('-s','--size',
		dest='size',
		default = "150,400",
		help="Filter primers by PRODUCT size.default [150,400].")

parser.add_option('-d','--dist',
		dest='dist',
		default = 4,
		help='Filter param of hairpin, which means distance of the minimal paired bases. Default: 4. Example:(number of X) AGCT[XXXX]AGCT')

parser.add_option('-k','--kmer',
		dest='kmer',
		default=7,
		type="int",
		help='Filter param of middle dimer. This param is used to cut primer into kmers and step is 1, length of kmer is [default: 7].\n \
		This param also control the end dimer length, max end dimer length = [kmer]-1.')

parser.add_option('-l','--loss',
                dest='loss',
                default = 4,
		type="int",
                help='Filter param of hairpin, which means distance of the minimal paired bases. Default: 4. Example:(number of X) AGCT[XXXX]AGCT')

parser.add_option('-a','--adaptor',
		dest='adaptor',
		default="0",
		type="int",
		help='Adaptor length. \n \
		>F \n \
		TCTTTCCCTACACGACGCTCTTCCGATCT \n \
		>R \n \
		TGGAGTTCAGACGTGTGCTCTTCCGATCT')

parser.add_option('-w','--window',
		dest='window',
		default = "10",
		help='Filter primers by primer distance. sliding window of primers. the minimal distance of primer_Fs or primer_Rs should greater than [window]')

parser.add_option('-o','--out',
		dest='out',
		help='Prefix of out file: candidate primers')

(options,args) = parser.parse_args()
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
count = 3   #times for the packages install
while count:
	try:
		import Bio   #
		print ('Dependent package Biopython is OK.\nDpendent module Bio is OK.')
		break
	except:
		print ('Dependent package Biopython is not found!!! \n Start intalling ....')
		os.system('pip install biopython')
		count -= 1
		continue

import collections
from collections import defaultdict
import time
from time import strftime
import re
import os
import os.path
import pandas as pd
import Bio
from Bio.Seq import Seq

matchnumber = int(options.num)
seqnumber = options.ref + ".number"
os.system("grep '>' {} | wc -l | cut -f 1 > {}".format(options.ref,seqnumber))
with open(seqnumber,"r") as f:
	seq_number = int(f.readline().strip())
print("INFO {} Start: Loading sequence......\n \nNumber of total sequence: {}.\n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time())),seq_number))
degeprimer = open(options.input,"r")
out = open(options.out,"w")
if seq_number > 10:
	frac = float(options.frac)
else:
	frac = float(options.frac) * 0.7
print("INFO {}: Load primers from DEGEPRIMER......\n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
########################## step1 screen #############################
#########################################################################################
# replace degenerate base to [A,C,G,T].
import math  # 
from functools import reduce
from operator import mul  # 

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
#################################################################
distance = int(options.dist)
#4bp revcomplement, distance = distane
def hairpin_check(primer):
	n = 0
	check = "FALSE"
	while n <= len(primer)-5-5-distance:
		kmer = dege_trans(primer[n:n+5])
		left = dege_trans(primer[n+5+distance:])
		for k in kmer:
			for l in left:
				if re.search(str(Seq(k).reverse_complement()),l):  #complement is ok, but reverse_complement must throw away.
					check = "TRUE"
					break
			if check == "TRUE":
				break
		if check == "TRUE":
			break
		n += 1

	if check == "TRUE":
		return True
	else:
		return False

###########################################################

pre_primer_pos = {}
pre_primer_frac = {}
pre_primer_match = {}
pre_primer_GC = {}
gc_content = options.gc.split(",")
print("GC content filter:{}".format(gc_content))
for i in degeprimer:
	if i.startswith("Pos"):
		pass
	else:
		i = i.strip().split("\t")
		position = int(i[0])
		number_match = int(i[6])
		fraction = float(number_match/seq_number)
		primer = i[5]
		primer_len = len(primer)
		GC_content = round((list(primer).count("G") + list(primer).count("C"))/len(list(primer)),3)
		if fraction < frac:
			pass
		elif re.search("AAAA|TTTT|CCCC|GGGG|AGAGAG|ATATAT|ACACAC|CTCTCT|CACACA|CGCGCG|TATATA|TCTCTC|TGTGTG|GAGAGA|GTGTGT|GCGCGC",primer):
			pass
		elif GC_content > float(gc_content[1]) or GC_content < float(gc_content[0]):
			pass
			#print("GC content: {}; pass".format(GC_content))
		elif number_match < matchnumber:
			pass
		elif hairpin_check(primer):
			pass
		else:
			#print(primer)
			pre_primer_match[primer] = number_match
			pre_primer_pos[primer] = position
			pre_primer_frac[primer] = fraction
			pre_primer_GC[primer] = GC_content
degeprimer.close()
### sort the primer dictionary and check the max fraction
#pre_primer_pos = sorted(pre_primer_pos.items(), key = lambda k:k[1], reverse = True)
#sorted(pre_primer_frac.items(), key = lambda k:k[1],reverse=True)
#print(pre_primer_frac)
### PRODUCT size >= 150
if bool(pre_primer_pos):
	maxpos = max(pre_primer_pos.values())
	minpos = min(pre_primer_pos.values())
	if maxpos - minpos < 150:
		print("\n****************\nthe smallest site is {} and the largest site is {} \nthe product size is smaller than 150 and exit\n****************\n".format(minpos,maxpos))
	else:
		print("Non-filter primers number:{}".format(len(pre_primer_frac)))
		for i in pre_primer_frac.keys():
			out.write(">"+str(pre_primer_pos[i]) + "\n" + i + "\n")
			#out.write(">"+str(pre_primer_pos[i]) + "|match: " + str(pre_primer_match[i]) + "|fraction: " + str(round(pre_primer_frac[i],2)) + "|GC_content: " + str(pre_primer_GC[i]) + "\n" + i + "\n")
	out.close()

else:
	print("Non-cancidate primer for the input sequence\n")
################################ blast ###################################
size = os.path.getsize(options.out)
try:
	1/size
except:
	print("Non-primers!") 
	raise SystemExit()	

if re.search("/", options.ref):
	db_dir = options.ref.split("/")
	db = options.ref + ".db/" + db_dir[-1]
else:
	db = options.ref + ".db/" + options.ref

if os.path.isdir(db):
	print("Databse for blast is OK!!!\n")
else:
	#os.system("mkdir -p {}".format(db))
	print("INFO {}: make blastdb......\n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
	os.system("makeblastdb -dbtype nucl -in {} -out {}".format(options.ref,db))

blast_out = options.ref + ".blastout"

print("INFO {}: blastn with 2 mismatich (degenerate);\n\nPrimer length is: {}.\nPerfect match len: {}.\n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time())),primer_len,primer_len-2))

if os.path.isfile(blast_out):
	print("\nBlast is done! Please check your blastout file, make sure blastout is ok!!!\n")
else:
	os.system("blastn -task megablast -word_size {} -outfmt '6 qseqid qlen qstart qend saccver slen sstart send bitscore length pident' -query {} -db {} -num_threads 20 -out {} -max_target_seqs 20000000".format(primer_len-2,options.out,db,blast_out))

######################### primerID : primerSeq ###########################
primer_seq = {}
with open(options.out,'r') as primers:
	for p in primers:
		if p.startswith(">"):
			key = p.lstrip(">").strip()
		else:
			primer_seq[int(key)] = p.strip()
###########################################################################
#dna_seq.complement()
#dna_seq.reverse_complement()
#########################  dimer check  ################################
kmer = int(options.kmer)
def Penalty_points(length, GC, d1, d2):
	return math.log((2 ** length * 2 ** GC) / ((d1 + 0.2) * (d2 + 0.1)), 10)
loss = options.loss
adaptor_len = options.adaptor
#primer_end_set = set()
def dimer_check(primer_F,primer_R): #Caution: primer_key_set must be a global var.
	check = "F"
	primer_F_R = set(dege_trans(primer_F)).union(set(dege_trans(primer_R)))
	small_end_set = set()
	n = 0
	while n < kmer - 4:
		four_F_seq = dege_trans(primer_F[-n-4:])
		four_R_seq = dege_trans(primer_R[-n-4:])
		small_end_set = small_end_set.union(set(four_F_seq))
		small_end_set = small_end_set.union(set(four_R_seq))
		n += 1
	for seq in small_end_set:
		for primer in primer_F_R:
			end_length = len(seq)
			end_GC = seq.count("G") + seq.count("C")
			end_d1 = 0
			if re.search(str(Seq(seq).reverse_complement()), primer) or str(Seq(seq).complement()) in primer:
				if re.search(str(Seq(seq).reverse_complement()), primer):
					r_end = min((len(primer) - len(seq) - primer.index(str(Seq(seq).reverse_complement()))) \
						,primer.index(str(Seq(seq).reverse_complement()))+adaptor_len)
				else:
					r_end = 10000
				if str(Seq(seq).complement()) in primer:
					c_end = min((len(primer) - len(seq) - primer.index(str(Seq(seq).complement()))) \
						,primer.index(str(Seq(seq).complement())))
				else:
					c_end = 10000
				end_d2 = min(r_end,c_end)
				Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
				if Loss > 2: #min loss: length=4, GC=2, end_d1=0, end_d2>=4
						#       length=5, GC=2, end_d1=0, end_d2>=7
					check = "T"
		if check == "T":
			break
	if check == "T":
		pass
	else:
		n,k = 0,0
		middle_seq = {}
		primer_F_len,primer_R_len = len(primer_F),len(primer_R)
		while n < primer_F_len - kmer: 
			for f in range(len(primer_F)):
				for i in dege_trans(primer_F[n:n + kmer]):
					l_middle_seq = primer_F_len - n - kmer + 1
					r_middle_seq = n + adaptor_len
					middle_seq[i]= min(l_middle_seq,r_middle_seq)
				f += 1
			n += 1
		while k < primer_R_len - kmer:
			for r in range(len(primer_R)):
				for j in dege_trans(primer_R[k:k + kmer]):
					l_middle_seq = primer_R_len - n - kmer + 1
					r_middle_seq = k + adaptor_len	
					middle_seq[j] = min(l_middle_seq,r_middle_seq)
				r += 1
			k += 1
		for seq in middle_seq.keys():
			for R in primer_F_R:
				if re.search(str(Seq(seq).reverse_complement()), R) or re.search(str(Seq(seq).complement()), R):
					middle_length = len(seq)
					middle_GC = seq.count("G") + seq.count("C")
					middle_d1 = middle_seq[seq]
					if re.search(str(Seq(seq).reverse_complement()), R):
						r_middle = min((len(R) - len(seq) - R.index(str(Seq(seq).reverse_complement()))),
							R.index(str(Seq(seq).reverse_complement()))+adaptor_len)
					else:
						r_middle = 10000
					if re.search(str(Seq(seq).complement()), R):
						c_middle = min((len(R) - len(seq) -R.index(str(Seq(seq).complement()))),  
							R.index(str(Seq(seq).complement()))+adaptor_len)
					else:
						c_middle = 10000
					middle_d2 = min(r_middle,c_middle)
					Loss = Penalty_points(middle_length, middle_GC, middle_d1, middle_d2)
					if Loss > loss: #min loss: length=7, GC=3, d1=2, d2=6
							#	   length=7, GC=3, d1=3, d2=4
						check = "T"
						break
			if check == "T":
				break
	if check == "T": #contain end dimer
		return True
	else:
		return False


#########################  paired primers  ################################
user_product_size = options.size.split(",")
print("PCR PRODUCT SIZE: {} - {}".format(user_product_size[0],user_product_size[1]))

window = int(options.window)
primer2speciesID_start_dict = defaultdict(list)
primer2specify_mismatch_dict = defaultdict(list)
def get_PCR_PRODUCT(blastout,output,candidate_primer_out,candidate_primer_txt):
	tb_out = pd.DataFrame(columns=["speciesID","primer_F:R_dege","primer_F:R_seq","pF_start","pR_end","product_size","pF_fraction","pR_fraction","pF_target_number","pR_target_number","f_mismatch","r_mismatch"])
	last_primer_start = 0
	tb_primer  = pd.DataFrame(columns=["primer_F:R_dege","primer_F:R_seq"])
	for idx,row in blastout.iterrows():
		primerID = int(row[0])
		specisesID = row[4]
		product_start = int(row[6])
		mismatch = len(primer_seq[primerID]) - int(row[9])
		value = (specisesID, product_start)
		#value: (NCBI_ID: mapped_position)
		primer2speciesID_start_dict[primerID].append(value)
		#primer position: (NCBI_ID: mapped_position)
		key = (primerID,specisesID)
		primer2specify_mismatch_dict[key].append(mismatch)
	print("primer information is OK !!!")
	for primer_start in primer2speciesID_start_dict.keys():
		if primer_start - last_primer_start <= window: #primer F dist >= 6
			pass
		else:
			last_primer_end = 0
			for primer_end in primer2speciesID_start_dict.keys():
				if primer_end - primer_start > int(user_product_size[1]): #prunding
					break
				elif primer_end - last_primer_end <= window:
					pass
				else:
					if primer_end - primer_start >= int(user_product_size[0]) and primer_end - primer_start <= int(user_product_size[1]):
						primer_F_R = str(primer_start) + ":" + str(primer_end)
						primer_F_seq = primer_seq[primer_start]
						primer_R_seq = str(Seq(primer_seq[primer_end]).reverse_complement())
						if dimer_check(primer_F_seq,primer_R_seq): #remove dimer
							print("Dimer primer: {} - {}. Removing...".format(primer_F_seq,primer_R_seq))
							pass
						else:
							primer_F_R_seq = primer_seq[primer_start] + ":" + str(Seq(primer_seq[primer_end]).reverse_complement())
							primer_len = len(primer_seq[primer_start])
							primer_F_dict = dict(primer2speciesID_start_dict[primer_start])
							primer_R_dict = dict(primer2speciesID_start_dict[primer_end])
							primer_F_target = len(primer_F_dict)
							primer_R_target = len(primer_R_dict)
							for ID in primer_F_dict.keys(): # NCBI_ID in primer_start
								if ID not in primer_R_dict.keys(): # NCBI_ID in primer_end
									pass
								else:	
									speciesID = ID
									primer_F_specify = len(primer2specify_mismatch_dict[(primer_start,speciesID)])
									primer_R_specify = len(primer2specify_mismatch_dict[(primer_end,speciesID)])
									#print(primer_F_specify)
									if primer_F_specify > 1 or primer_R_specify > 1:
										pass
									else:
										product_start = int(primer_F_dict[ID])
										product_stop = int(primer_R_dict[ID]) + primer_len
										product_size = int(primer_R_dict[ID]) + primer_len - int(primer_F_dict[ID]) + 1
										primer_F_mismatch = primer2specify_mismatch_dict[(primer_start,speciesID)]
										primer_R_mismatch = primer2specify_mismatch_dict[(primer_end,speciesID)]
										primer_F_fraction = pre_primer_frac[primer_seq[primer_start]]
										primer_R_fraction = pre_primer_frac[primer_seq[primer_end]]
										tb_out = tb_out.append({
											"speciesID":speciesID,
											"primer_F:R_dege":primer_F_R,
											"primer_F:R_seq":primer_F_R_seq,
											"pF_start":product_start,
											"pR_end":product_stop,
											"product_size":product_size,
											"pF_fraction":primer_F_fraction,
											"pR_fraction":primer_R_fraction,
											"pF_target_number":primer_F_target,
											"pR_target_number":primer_R_target,
											"f_mismatch":primer_F_mismatch,
											"r_mismatch":primer_R_mismatch
										},ignore_index=True)
					last_primer_end = primer_end
					last_primer_start = primer_start	
	tb_out.sort_values(by=["pF_target_number","pR_target_number"],inplace=True,ascending=False)
	tb_out.to_csv(output,index=False,sep="\t")
	primer_txt = tb_out.loc[:,["primer_F:R_dege","primer_F:R_seq","pF_target_number","pR_target_number"]].drop_duplicates()
	candidate_primer_txt.write(options.out)
	for idx,row in primer_txt.iterrows():
		product_start_end = row[0].split(":")
		product_len = str(int(product_start_end[1]) -int(product_start_end[0]))
		primer_sequence = row[1].split(":")
		targets_number = str(min(row[2],row[3]))
		#file for hairpin re-check
		candidate_primer_out.write(">" + product_start_end[0] + "_F" +"\n" + primer_sequence[0]+"\n>" + product_start_end[1] + "_R" +"\n" + primer_sequence[1]+"\n")
		#file for multiPCR primer select
		candidate_primer_txt.write("\t" +primer_sequence[0]+"\t" + primer_sequence[1] + "\t" + product_len + "\t" + targets_number + "\t" + row[0])
	candidate_primer_txt.write("\n")
#primer_check
paired = options.out.rstrip(".candidate.primers.txt") + ".paired.nt.Check"
paired_primer = open(paired,"w")
#blastn output
blast_results = open(blast_out,"r")
results = pd.read_table(blast_results)
candidate_primer = options.out.rstrip(".candidate.primers.txt") + ".candidate.primers.fa"
candidate_primer_for_set = options.out
candidate_primer_fa = open(candidate_primer,"w")
candidate_primer_txt = open(candidate_primer_for_set,"w")
get_PCR_PRODUCT(results,paired_primer,candidate_primer_fa,candidate_primer_txt)
blast_results.close()
paired_primer.close()
candidate_primer_fa.close()


print("INFO {}: Done!!!\n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
#########################  Done!  ################################

