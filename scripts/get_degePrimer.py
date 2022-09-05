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
			Options: {-f -n [num] -s [250,500] -g [0.4,0.6] -d [ditance] -a ","}.')
parser.add_option('-i', '--input',
                  dest='input',
                  help='Input file: degeprimer out.')

parser.add_option('-r', '--ref',
                  dest='ref',
                  help='Reference sequence file: all the sequence in 1 fasta, for example: (Cluster_96_171.fa).')

parser.add_option('-n', '--num',
                  dest='num',
                  default="500",
                  type="int",
                  help='Filter primers by match rank: sort candidate primers by match number and select top {N} for '
                       'next steps. Default: 500.')

parser.add_option('-g', '--gc',
                  dest='gc',
                  default="0.4,0.65",
                  help="Filter primers by GC content. Default [0.4,0.6].")

parser.add_option('-f', '--fraction',
                  dest='fraction',
                  default="0.6",
                  type="float",
                  help="Filter primers by match fraction. Default: 0.6.")

parser.add_option('-t', '--term',
                  dest='term',
                  default="3",
                  type="int",
                  help="Filter primers by degenerate base position. e.g. [-t 4] means I dont want degenerate base "
                       "appear at the term four base when primer pre-filter. Default: 4.")

parser.add_option('-p', '--position',
                  dest='position',
                  default="9",
                  type="int",
                  help="Filter primers by mismatch position. e.g. [-p 8] means I dont want mismatch appear  at the "
                       "term eight base when primer checking. Default: 8.")

parser.add_option('-s', '--size',
                  dest='size',
                  default="250,500",
                  help="Filter primers by PRODUCT size. Default [250,500].")

parser.add_option('-d', '--dist',
                  dest='dist',
                  default=4,
                  help='Filter param of hairpin, which means distance of the minimal paired bases. Default: 4. '
                       'Example:(number of X) AGCT[XXXX]AGCT.')

parser.add_option('-a', '--adaptor',
                  dest='adaptor',
                  default="TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT",
                  type="str",
                  help='Adaptor sequence, which is used for NGS next. Hairpin or dimer detection for adaptor--primer. '
                       '\n \ For example: TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT (Default). If '
                       'you dont want adaptor, use [","] ')

parser.add_option('-m', '--maxseq',
                  dest='maxseq',
                  default=500,
                  type="int",
                  help='Limit of sequence number. Default: 500')

parser.add_option('-o', '--out',
                  dest='out',
                  help='Prefix of out file: candidate primers')

(options, args) = parser.parse_args()

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
from statistics import mean
# import operator

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
count = 3  # times for the packages install
while count:
    try:
        import Bio  #

        print('Dependent package Biopython is OK.\nDpendent module Bio is OK.')
        break
    except:
        print('Dependent package Biopython is not found!!! \n Start intalling ....')
        os.system('pip install biopython')
        count -= 1
        continue

########################## step1 screen #############################
#########################################################################################
# replace degenerate base to [A,C,G,T].
import math  # 
from functools import reduce
from operator import mul  # 

'''
def get_dict_key(d_dict, value):
	d_key = list(d_dict.keys())[list(d_dict.values()).index(value)]  # not for multi keys corresponding 1 value
	return d_key
'''
degenerate_pair = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

degenerate_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
                    "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}


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

# 5 bp revcomplement, distance = distane
def hairpin_check(primer):
    n = 0
    check = "FALSE"
    while n <= len(primer) - 5 - 5 - distance:
        kmer = dege_trans(primer[n:n + 5])
        left = dege_trans(primer[n + 5 + distance:])
        for k in kmer:
            for l in left:
                if re.search(str(Seq(k).reverse_complement()),
                             l):  # complement is ok, but reverse_complement must throw away.
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

def GC_fraction(sequence):
    sequence_expand = dege_trans(sequence)
    GC_list = []
    for seq in sequence_expand:
        GC_list.append(round((list(seq).count("G") + list(seq).count("C")) / len(list(seq)), 3))
    GC_average = mean(GC_list)
    return GC_average


###########################################################

def pre_filter(degeprimer, GC, maxseq, frac, rank_number,tmp):
    # primers pre-filter.
    pre_primer_pos = {}
    global pre_primer_frac
    pre_primer_match = {}
    pre_primer_GC = {}
    gc_content = GC
    with open(degeprimer, "r") as degeprimer:
        degeprimer_results = pd.read_table(degeprimer, header="infer")
        degeprimer_results.sort_values(by=["NumberMatching"], inplace=True, ascending=False)
        # loc []; iloc [). Suggestion: label selection use loc, index selection use iloc.
        degeprimer_filter = pd.DataFrame(degeprimer_results.iloc[0:rank_number])
        degeprimer_filter.sort_values(by=["Pos"], inplace=True, ascending=True)
        print(degeprimer_filter)
        # print(degeprimer_filter.shape)
        for idx, row in degeprimer_filter.iterrows():
            position = int(row[0])
            number_match = int(row[6])
            if seq_number <= maxseq:
                fraction = float(number_match / seq_number)
            else:
                fraction = float(number_match / maxseq)
            primer = row[5]
            GC_content = GC_fraction(primer)
            if fraction < frac:
                pass
            elif re.search("AAAA|CCCC|GGGG|TTTT", primer):
                pass
            # elif re.search("CCC|CCG|CGG|CGC|GCC|GCG|GGC|GGG",primer[-3:]):
            # pass
            elif GC_content > float(gc_content[1]) or GC_content < float(gc_content[0]):
                pass
            else:
                pre_primer_match[primer] = number_match
                pre_primer_pos[primer] = position
                pre_primer_frac[primer] = fraction
                pre_primer_GC[primer] = GC_content

    ##### check PRODUCT size #####
    if bool(pre_primer_pos):
        maxpos = max(pre_primer_pos.values())
        minpos = min(pre_primer_pos.values())
        if maxpos - minpos < int(user_product_size[0]):
            print("\n***********\nthe smallest site is {} and the largest site is {}. \n \
             Product size is smaller than min product size and exit\n***********\n".format(minpos, maxpos,
                                                                                           user_product_size[0]))
            raise SystemExit()
        else:
            print("\n Non-filter primers number: {}. \n".format(len(pre_primer_frac)))
            with open(tmp, "w") as out:
                for i in pre_primer_frac.keys():
                    out.write(">" + str(pre_primer_pos[i]) + "\n" + i + "\n")
    else:
        print("Non-cancidate primer for the input sequence.\n")
        raise SystemExit()
    ###### check primer #####
    size = os.path.getsize(tmp)
    try:
        1 / size
    except:
        print("Non-primers!")
        raise SystemExit()


######################### MAP: bowtie2 ###########################
def map(path_2_ref, tmp_expand, sam_out, sam_for_out):
    # path_2_ref is used to build bowtie2 index
    # path_2_out is used to output samfile
    title = path_2_ref.split("/")
    db = path_2_ref + ".db/"
    if os.path.isdir(db):
        print("Databse for mapping is OK!!!\n")
    else:
        os.system("mkdir -p {}".format(db))
        print("INFO {}: Indexing......\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))
        os.system("bowtie2-build {} {}/{}".format(path_2_ref, db, title[-1]))
    db_title = db + "/" + title[-1]
    if os.path.isfile(sam_out):
        print("Mapping is done! Please check your map file, make sure sam_out is ok!!!\n")
    else:
        os.system(
            # -L  length of seed substrings
            "bowtie2 -N 1 -L 31 -a -x {} -f -U {} -S {}".format(
                db_title, tmp_expand, sam_out))
        os.system(
            # remove primer which matched to the reverse strand
            "samtools view -F 16 {} > {}".format(
                sam_out, sam_for_out))


######################### primerID : primerSeq ###########################
def primerSeq(pre_primer_seq):
    primer_seq = {}
    with open(pre_primer_seq, 'r') as primers:
        for p in primers:
            if p.startswith(">"):
                key = p.lstrip(">").strip()
            else:
                primer_seq[int(key)] = p.strip()
    return primer_seq


#########################  dimer check  ################################

def Penalty_points(length, GC, d1, d2):
    return math.log((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)), 10)


def current_end(primer_F, primer_R):
    primer_F_extend = adaptor[0] + primer_F
    primer_R_extend = adaptor[1] + primer_R
    primer_F_len = len(primer_F)
    primer_R_len = len(primer_R)
    end_seq = set()
    for a in range(primer_F_len - 5):
        F_end_seq = dege_trans(primer_F_extend[-a - 5:])
        end_seq = end_seq.union(set(F_end_seq))
        F_start_seq = dege_trans(primer_F_extend[:a + 5])
        end_seq = end_seq.union(set(F_start_seq))
        a += 1
    for b in range(primer_R_len - 5):
        R_end_seq = dege_trans(primer_R_extend[-b - 5:])
        end_seq = end_seq.union(set(R_end_seq))
        R_end_seq = dege_trans(primer_R_extend[:b + 5])
        end_seq = end_seq.union(set(R_end_seq))
        b += 1
    return end_seq


def dimer_check(primer_F, primer_R):  # Caution: primer_key_set must be a global var.
    check = "F"
    primer_F_R = set(dege_trans(primer_F)).union(set(dege_trans(primer_R)))
    small_end_set = current_end(primer_F, primer_R)
    for seq in small_end_set:
        for primer in primer_F_R:
            end_length = len(seq)
            end_GC = seq.count("G") + seq.count("C")
            end_d1 = 0
            if re.search(str(Seq(seq).reverse_complement()), primer):
                end_d2 = min((len(primer) - len(seq) - primer.index(str(Seq(seq).reverse_complement())))
                             , primer.index(str(Seq(seq).reverse_complement())) + adaptor_len)
                Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
            else:
                end_d2 = "NA"
                Loss = 0
            if Loss > 3:
                check = "T"
                print("Dimer seq: {}; dimer length: {}; dimer position: {}; primer: {}".format(
                    str(Seq(seq).reverse_complement()),
                    str(len(str(Seq(seq).reverse_complement()))), str(end_d2), primer))
                break
        if check == "T":
            break
    if check == "T":
        return True
    else:
        return False


#########################  paired primers  ################################

def get_PCR_PRODUCT(sam, output, candidate_primer_out, candidate_primer_txt,dege_pos,bed_file):
    primer2speciesID_start_dict = defaultdict(list)
    primer2specify_mismatch_dict = {}
    tb_out = pd.DataFrame(columns=["speciesID", "primer_F:R_dege", "primer_F:R_seq", "pF_start", "pR_end",
                                   "product_size", "pF_fraction", "pR_fraction", "pF_target_number",
                                   "pR_target_number", "f_mismatch_position", "r_mismatch_position"])
    #### primer information ####
    results = pd.read_table(sam, header=None)
    # print(results.iloc[0])
    for rank in range(11,len(results.iloc[0])):
        # print(results.iloc[0][rank])
        if re.search("MD:Z:", results.iloc[0][rank]):
            p_rank = rank
    # mismatch_pattern = re.compile('XM:i:(\d+)')
    position_pattern = re.compile("[A-Z]?(\d+)")
    for idx, row in results.iterrows():
        primerID = int(row[0].split("_")[0])
        specisesID = row[2]
        product_start = int(row[3])
        # XM:i mismatch number
        # mismatch = mismatch_pattern.search(row[14]).group(1)
        position = position_pattern.search(row[p_rank][-2:]).group(1)
        # distance between mismatch position and 3' term must greater than 8
        if int(position) < mismatch_position:
            pass
        else:
            mis_key = (primerID, specisesID)
            value = (specisesID, product_start)
            if mis_key in primer2specify_mismatch_dict.keys():
                # primer have more than one targets sites. The match is equal to or higher than any other match.
                if position > primer2specify_mismatch_dict[mis_key]:
                    primer2specify_mismatch_dict[mis_key] = position
                    primer2speciesID_start_dict[primerID].append(value)
                else:
                    primer2speciesID_start_dict[primerID].insert(0,value)
            else:
                primer2speciesID_start_dict[primerID].append(value)
                primer2specify_mismatch_dict[mis_key] = position
    primer2speciesID_start_dict = collections.OrderedDict(primer2speciesID_start_dict)
    print("primer information is OK !!!")
    #### paired primers ####
    for primer_start in primer2speciesID_start_dict.keys():
        position = list(primer2speciesID_start_dict.keys()).index(primer_start)
        for primer_end in list(primer2speciesID_start_dict.keys())[position:]:
            if primer_end - primer_start + len(primer_seq[primer_end]) < int(user_product_size[0]):
                pass
            elif primer_end - primer_start + len(primer_seq[primer_end]) > int(user_product_size[1]):
                break
            else:
                primer_F_R = str(primer_start) + ":" + str(primer_end + len(primer_seq[primer_end]))
                primer_F_seq = primer_seq[primer_start]
                primer_R_seq = str(Seq(primer_seq[primer_end]).reverse_complement())
                if dege_filter_in_term_N_bp(primer_F_seq,dege_pos) or dege_filter_in_term_N_bp(primer_R_seq,dege_pos):
                    print("Degenerate base detected in the term {} nucleotides. removing...".format(dege_pos))
                elif dimer_check(primer_F_seq, primer_R_seq):  # remove dimer
                    print("Dimer primer detected: {} - {}. Removing...".format(primer_F_seq, primer_R_seq))
                elif hairpin_check(adaptor[0] + primer_F_seq) or hairpin_check(adaptor[1] + primer_R_seq):
                    print("Hairpin structure detected. removing...")
                else:
                    primer_F_R_seq = primer_seq[primer_start] + ":" + str(
                        Seq(primer_seq[primer_end]).reverse_complement())
                    primer_F_dict = dict(primer2speciesID_start_dict[primer_start])
                    primer_R_dict = dict(primer2speciesID_start_dict[primer_end])
                    for ID in primer_F_dict.keys():  # NCBI_ID in primer_start
                        if ID not in primer_R_dict.keys():  # NCBI_ID in primer_end
                            pass
                        else:
                            speciesID = ID
                            product_start = int(primer_F_dict[ID])-1
                            product_stop = int(primer_R_dict[ID]) + len(primer_seq[primer_end])-1
                            product_size = int(primer_R_dict[ID]) + len(primer_seq[primer_end]) - int(
                                primer_F_dict[ID]) + 1
                            primer_F_mismatch = primer2specify_mismatch_dict[(primer_start, speciesID)]
                            primer_R_mismatch = primer2specify_mismatch_dict[(primer_end, speciesID)]
                            primer_F_fraction = round(len(primer_F_dict) / seq_number, 2)
                            primer_R_fraction = round(len(primer_R_dict) / seq_number, 2)
                            primer_F_target = len(primer_F_dict)
                            primer_R_target = len(primer_R_dict)
                            bed_file.write(speciesID + "\t" + str(product_start) + "\t" + str(product_stop) +"\n")
                            tb_out_local = pd.DataFrame({
                                "speciesID": speciesID,
                                "primer_F:R_dege": primer_F_R,
                                "primer_F:R_seq": primer_F_R_seq,
                                "pF_start": product_start,
                                "pR_end": product_stop,
                                "product_size": product_size,
                                "pF_fraction": primer_F_fraction,
                                "pR_fraction": primer_R_fraction,
                                "pF_target_number": primer_F_target,
                                "pR_target_number": primer_R_target,
                                "f_mismatch_position": primer_F_mismatch,
                                "r_mismatch_position": primer_R_mismatch
                            }, index=[0])
                            tb_out = pd.concat([tb_out, tb_out_local], axis=0, ignore_index=True)
    tb_out.sort_values(by=["pF_target_number", "pR_target_number"], inplace=True, ascending=[False, False])
    # print(tb_out)
    tb_out.to_csv(output, index=False, sep="\t")
    primer_txt = tb_out.loc[:,["primer_F:R_dege", "primer_F:R_seq", "pF_target_number", "pR_target_number"]].drop_duplicates()
    primer_txt.sort_values(by=["pF_target_number", "pR_target_number"],inplace=True,ascending=[False, False])
    # print(primer_txt)
    candidate_primer_txt.write(options.out)
    for idx, row in primer_txt.iterrows():
        product_start_end = row[0].split(":")
        product_len = str(int(product_start_end[1]) - int(product_start_end[0]))
        primer_sequence = row[1].split(":")
        targets_number = str(min(row[2], row[3]))
        # file for hairpin re-check
        candidate_primer_out.write(
            ">" + product_start_end[0] + "_F" + "\n" + primer_sequence[0] + "\n>" + product_start_end[1] + "_R" + "\n" +
            primer_sequence[1] + "\n")
        # file for multiPCR primer select
        candidate_primer_txt.write(
            "\t" + primer_sequence[0] + "\t" + primer_sequence[1] + "\t" + product_len + "\t" + targets_number + "\t" +
            row[0])
    candidate_primer_txt.write("\n")

#########################  degenerate primers  ################################

def fa_dege_trans(fa, output):
    out = open(output, "w")
    with open(fa, "r") as dege_fa:
        for i in dege_fa:
            if i.startswith(">"):
                Id = i.strip()
            else:
                sequence = i.strip()
                expand_seq = dege_trans(sequence)
                if len(expand_seq) > 1:
                    for j in range(len(expand_seq)):
                        ID = Id + "_" + str(j)
                        out.write(ID + "\n" + expand_seq[j] + "\n")
                else:
                    ID = Id + "_0"
                    out.write(ID + "\n" + expand_seq[0] + "\n")
    out.close()

#########################  degenerate position filter  ################################

def dege_filter_in_term_N_bp(sequence, term):
    if term == 0:
        term_base = ["A"]
    else:
        term_base = sequence[-term:]
    score = score_trans(term_base)
    if score > 1:
        return True
    else:
        return False

#########################  main  ################################

if __name__ == "__main__":
    #### hairpin distance ####
    distance = int(options.dist)
    print("Hairpin distance: {}. \n".format(distance))
    #### adaptor ####
    adaptor = options.adaptor.split(",")
    adaptor_len = len(adaptor[0])
    print("Adaptor: {} - {}.\n".format(adaptor[0], adaptor[1]))

    #### Product size ####
    user_product_size = options.size.split(",")
    print("Product size: {} - {}.\n".format(user_product_size[0], user_product_size[1]))

    #### max_Seq ####
    max_seq = options.maxseq
    print("Max sequence number: {}.\n".format(max_seq))

    #### degenerate position ####
    dege_pos = options.term
    print("Users dont want degenerate base appear at: {}.\n".format(dege_pos))

    #### GC content ####
    GC = options.gc.split(",")
    print("GC content (candidate primer): {} - {}.\n".format(GC[0], GC[1]))

    #### rank_number ####
    rank_number = options.num
    print("Users choose top {} primers for the next steps.\n".format(rank_number))

    #### mismatch position ####
    mismatch_position = options.position
    print("Mismatch position should not located at term {} bp when primer checking.\n".format(mismatch_position))

    #### stastic of seq_number ####
    sequence_number = options.ref + ".number"
    os.system("grep '>' {} | wc -l | cut -f 1 > {}".format(options.ref, sequence_number))
    with open(sequence_number, "r") as f:
        seq_number = int(f.readline().strip())
    print("INFO {} Start: Loading sequence......\n Number of total sequence: {}.\n".format(
        time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())), seq_number))
    if seq_number > 10:
        frac = float(options.fraction)
    else:
        frac = float(options.fraction) * 0.7

    #### pre_filter ####
    pre_primer_frac = {}
    tmp = options.out.rstrip(".candidate.primers.txt") + ".pre_filter_primers.fa"
    pre_filter(options.input, GC, max_seq, frac, rank_number,tmp)

    #### tmp file ####
    tmp_expand = options.out.rstrip(".candidate.primers.txt") + ".pre_filter_primers.expand.fa"
    fa_dege_trans(tmp, tmp_expand)
    primer_seq = primerSeq(tmp)

    #### mapping with bowtie2 ####
    sam_out = options.out + ".sam"
    sam_for_out = options.out + ".for.sam"
    map(options.ref, tmp_expand, sam_out, sam_for_out)

    # primer_check
    paired = options.out.rstrip(".candidate.primers.txt") + ".paired.Check"
    paired_primer = open(paired, "w")
    sam_results = open(sam_for_out, "r")
    candidate_primer = options.out.rstrip(".candidate.primers.txt") + ".candidate.primers.fa"
    candidate_primer_fa = open(candidate_primer, "w")
    candidate_primer_txt = open(options.out, "w")
    bed_file = options.out.rstrip(".candidate.primers.txt") + ".bed"
    bed = open(bed_file,"w")
    get_PCR_PRODUCT(sam_results, paired_primer, candidate_primer_fa, candidate_primer_txt,dege_pos,bed)
    bed.close()
    sam_results.close()
    paired_primer.close()
    candidate_primer_fa.close()
    candidate_primer_txt.close()
    os.system("rm {} {}".format(sam_out,sam_for_out))
    print("INFO {}: Done!!!\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))

#########################  Done!  ################################
