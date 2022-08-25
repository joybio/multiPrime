#!/bin/python
"""
Get the final primerset.
"""
__date__ = "2022-7-6"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import re
import os
import sys
from sys import argv
import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog -i [input] -o [output] \n \
                        Options: {-s [step] -m [T] -k [kmer]}',version = "%prog 0.0.4")
parser.add_option('-i', '--input',
                  dest='input',
                  help='Input file: primers')

parser.add_option('-a','--adaptor',
                dest='adaptor',
                default="/share/data3/yangjunbo/multiPCR/202206/20220620_primerset/Primers_set/adapter.fa",
                help='Apdaptor file (useless for now). Default:/share/data3/yangjunbo/multiPCR/202206/20220620_primerset/Primers_set/adapter.fa. \n \
                        >F \n \
                        TCTTTCCCTACACGACGCTCTTCCGATCT \n \
                        >R \n \
                        TGGAGTTCAGACGTGTGCTCTTCCGATCT')

parser.add_option('-s', '--step',
                  dest='step',
                  default=4,
                  type="int",
                  help='distance between primers; column number of primer1_F to primer2_F.')

parser.add_option('-l', '--loss',
                  dest='loss',
                  default=2,
                  type="int",
                  help='loss function. Filter dimer of primers')

parser.add_option('-m', '--maxi',
                  dest='maxi',
                  default="T",
                  type="str",
                  help='which method: maximal or maximum. If -m [T] use maximal; else maximum')

parser.add_option('-d', '--dist',
                  dest='dist',
                  default=5,
                  type="int",
                  help='Filter param of hairpin.')

parser.add_option('-k', '--kmer',
                  dest='kmer',
                  default=7,
                  type="int",
                  help='Filter param of middle dimer. This param is used to cut primer into kmers and step is 1, length of kmer is [default: 7].\n \
                This param also control the end dimer length, max end dimer length = [kmer]-1.')

parser.add_option('-o', '--out',
                  dest='out',
                  help='Prefix of out file: candidate primers')
(options, args) = parser.parse_args()

import re
import math
from operator import mul
from functools import reduce

kmer = 9
step = 4
distance = 5
count = 3  # times for the packages install
while count:
    try:
        import Bio  #
        print('Dependent package Biopython is OK.\nDpendent module Bio is OK.')
        break
    except:
        print('Dependent package Biopython is not found!!! \n Start intalling ....')
        os.system('pip install Bio')
        count -= 1
        continue
from Bio.Seq import Seq

maximal_check = False
primer_end_set = set()
primer_middle_set = set()
primer_set = set() 

current_inner_end = set()
current_inner_middle = set()

middle_seq = {} # this is a dict of middle seq and its distance to the primer end. only used the nearest distance.

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


# 5 bp reverse complement, distance = distane. 
# Its already done in the get_degePrimer.py. So we did not use this fucntion
def hairpin_check(primer):
    n = 0
    check = "FALSE"
    while n <= len(primer) - 10 - distance:
        fragment = dege_trans(primer[n:n + 5])
        left = dege_trans(primer[n + 5 + distance::])
        for k in fragment:
            for le in left:
                if re.search(str(Seq(k).reverse_complement()), le):
                    # complement is ok, but reverse_complement must throw away.
                    check = "TRUE"
                    # print(n, str(Seq(k).reverse_complement()), le)
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


def Penalty_points(length, GC, d1, d2):
    return math.log((2 ** length * 2 ** GC) / ((d1 + 1) * (d2 + 1)), 10)

loss = options.loss
def inner_dimer_check(primer_F, primer_R):
    m, n, k = 0, 0, 0
    check = "F"
    Loss = 0
    primer_F_len = len(primer_F)
    primer_R_len = len(primer_R)
    global current_inner_end, current_inner_middle, middle_seq
    end_seq = set()
    for a in range(kmer - 4):
        F_end_seq_forward = dege_trans(primer_F[-a - 4:])
        end_seq = end_seq.union(set(F_end_seq_forward))
#        F_end_seq_reverse = dege_trans(primer_F[:a + 4])
#        end_seq = end_seq.union(set(F_end_seq_reverse))
        a += 1
    for b in range(kmer-4):
        R_end_seq_forward = dege_trans(primer_R[-b - 4:])
        end_seq = end_seq.union(set(R_end_seq_forward))
#        R_end_seq_reverse =dege_trans(primer_R[b + 4:])
#        end_seq = end_seq.union(set(R_end_seq_reverse))
        b += 1
    primer_F = dege_trans(primer_F)
    primer_R = dege_trans(primer_R)
    primer_F_R = set(primer_F + primer_R)
    for seq in end_seq:
        for R in primer_F_R:
            if re.search(str(Seq(seq).reverse_complement()), R):  # or re.search(str(Seq(seq).complement()), R):
                end_length = len(seq)
                end_GC = seq.count("G") + seq.count("C")
                end_d1 = 0
                end_d2 = min((len(R) - len(seq) - R.index(str(Seq(seq).reverse_complement()))) \
                         ,R.index(str(Seq(seq).reverse_complement())))
                Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                if Loss > loss:
                    check = "T"
                    break
        if Loss > loss:
            check = "T"
            break
    if check == "T":  
        pass
    else:
        while n <= primer_F_len - kmer:
            for f in range(len(primer_F)):
                for i in dege_trans(primer_F[f][n:n + kmer]):
                    if i in middle_seq.keys() and primer_F_len - n - kmer + 1 > middle_seq[i]:
                        pass
                    else:
                        middle_seq[i] = primer_F_len - n - kmer + 1  
                f += 1
            n += 1
        while k <= primer_R_len - kmer:
            for r in range(len(primer_R)):
                for j in dege_trans(primer_R[r][n:n + kmer]):
                    if j in middle_seq.keys() and primer_R_len - n - kmer + 1 > middle_seq[j]:
                        pass
                    else:
                        middle_seq[j] = primer_R_len - n - kmer + 1
                r += 1
            k += 1
        for seq in middle_seq.keys():
            for R in primer_F_R:
                if re.search(str(Seq(seq).reverse_complement()), R):#  or re.search(str(Seq(seq).complement()), R):
                    middle_length = len(seq)
                    middle_GC = seq.count("G") + seq.count("C")
                    middle_d1 = middle_seq[seq]
                    middle_d2 = min((len(R) - len(seq) - R.index(str(Seq(seq).reverse_complement()))), \
                                R.index(str(Seq(seq).reverse_complement())))
                    Loss = Penalty_points(middle_length, middle_GC, middle_d1, middle_d2)
                    if Loss > loss:
                        check = "T"
                        break
            if Loss > loss:
                check = "T"
                break

    if check == "T":
        return True
    else:
        current_inner_end = end_seq
        current_inner_middle = middle_seq
        return False

def outer_dimer_check(primer_F, primer_R):
    check = "F"
    global primer_set
    ###############################################################
    for end in current_inner_end:
        for primer in primer_set:
            if re.search(str(Seq(end).reverse_complement()), primer): #or re.search(str(Seq(end).complement()), primer):
                end_length = len(end)
                end_GC = end.count("G") + end.count("C")
                end_d1 = 0
                if re.search(str(Seq(end).reverse_complement()), primer):
                    end_d2 = min((len(primer) - len(end) - primer.index(str(Seq(end).reverse_complement()))) \
                             ,primer.index(str(Seq(end).reverse_complement())))
                    #if re.search(str(Seq(end).complement()), primer):
                    #    tmp = len(primer) - len(end) - primer.index(str(Seq(end).complement()))
                    #    end_d2 = min(end_d2,tmp)
                #else:
                    #end_d2 = len(primer) - len(end) - primer.index(str(Seq(end).complement()))
                Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                if Loss > loss:
                    check = "T"
                    break
        if check == "T":
            break
    ###############################################################
    if check == "T":  #
        pass
    else:
        primer_F = dege_trans(primer_F)
        primer_R = dege_trans(primer_R)
        primer_F_R = set(primer_F + primer_R)
        for end in primer_end_set:
            for primer in primer_F_R:
                if re.search(str(Seq(end).reverse_complement()),
                             primer):# or re.search(str(Seq(end).complement()), primer):
                    end_length = len(end)
                    end_GC = end.count("G") + end.count("C")
                    end_d1 = 0
                    if re.search(str(Seq(end).reverse_complement()),primer):
                        end_d2 = min((len(primer) - len(end) - primer.index(str(Seq(end).reverse_complement()))) \
                                 ,primer.index(str(Seq(end).reverse_complement())))
                        #if re.search(str(Seq(end).complement()), primer):
                        #    tmp = len(primer) - len(end) - primer.index(str(Seq(end).complement()))
                        #    end_d2 = min(end_d2,tmp)
                    #else:
                        #end_d2 = len(primer) - len(end) - primer.index(str(Seq(end).complement())) 
                    Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                    if Loss > loss:
                        check = "T"
                        break
            if check == "T":
                break
    ###############################################################
    if check == "T":  #
        pass
    else:
        for middle in primer_middle_set:
            for primer in primer_set:
                if re.search(str(Seq(middle).reverse_complement()), primer):
                    middle_length = len(middle)
                    middle_GC = middle.count("G") + middle.count("C")
                    middle_d1 = middle_seq[middle]
                    middle_d2 = min((len(primer) - len(middle) - primer.index(str(Seq(middle).reverse_complement())))
                                ,primer.index(str(Seq(middle).reverse_complement())))
                    Loss = Penalty_points(middle_length, middle_GC, middle_d1, middle_d2)
                    if Loss > loss: # 9 bp, 2**9*2**4/((2+1)*(2+1))=16224/9;
                        #print(str(Seq(middle).reverse_complement()), primer, Loss)
                        #print(Loss)
                        check = "T"
                        break
            if check == "T":
                break
    ###############################################################
    if check == "T":
        pass
    else:
        for middle in current_inner_middle:
            for primer in primer_set:
                if re.search(str(Seq(middle).reverse_complement()), primer):
                    middle_length = len(middle)
                    middle_GC = middle.count("G") + middle.count("C")
                    middle_d1 = middle_seq[middle]
                    middle_d2 = min((len(primer) - len(middle) - primer.index(str(Seq(middle).reverse_complement())))
                                ,primer.index(str(Seq(middle).reverse_complement())))
                    Loss = Penalty_points(middle_length, middle_GC, middle_d1, middle_d2)
                    if Loss > loss:
                        check = "T"
                        #print(Loss,middle_length, middle_GC, middle_d1, middle_d2)
                        break
            if check == "T":
                break
    ###############################################################
    if check == "T":
        return True
    else:
        return False
###############################################################

def greedy_primers(primers, row_num, clique):  
    global primer_end_set, maximal_check, primer_middle_set, primer_set
    primer_end_set = set()
    primer_middle_set = set()
    jdict = {}
    primer_end_dict = {}  
    primer_middle_dict = {}  
    primer_dict = {}
    row_pointer = 0
    blank_row = 0
    while row_pointer < row_num:
        column_pointer = 1
        if maximal_check:
            break
        else:
            if len(primers[row_pointer]) <= 1:
                row_pointer += 1
                blank_row += 1
                # print(row_pointer)
            else:
                while column_pointer <= len(primers[row_pointer]) - step:
                    if maximal_check:
                        break
                    else:
                        if inner_dimer_check(primers[row_pointer][column_pointer],
                                             primers[row_pointer][column_pointer + 1]) or \
                                outer_dimer_check(primers[row_pointer][column_pointer],
                                                  primers[row_pointer][column_pointer + 1]):
                            # print(row_pointer, primers[row_pointer][column_pointer],
                            #      primers[row_pointer][column_pointer + 1])
                            column_pointer += step
                            while column_pointer > len(primers[row_pointer]) - step:
                                row_pointer -= 1
                                if row_pointer < blank_row:
                                    maximal_check = True
                                    print("Non maximum primer set. Try maximal primer set!")
                                    break
                                else:
                                    # print(row_pointer, column_pointer)
                                    column_pointer = jdict[row_pointer] + step
                                    primer_end_set = primer_end_dict[row_pointer]
                                    primer_middle_set = primer_middle_dict[row_pointer]
                                    primer_set = primer_dict[row_pointer]
                        else:
                            if column_pointer <= len(primers[row_pointer]) - step:
                                # print(row_pointer, column_pointer, len(primer_middle_set), primer_middle_set)
                                clique.append(
                                    primers[row_pointer][0] + "\t" + str(column_pointer) + "\t" +
                                        primers[row_pointer][column_pointer] + "\t" +
                                        primers[row_pointer][column_pointer + 1] + "\t" +
                                        primers[row_pointer][column_pointer + 2] + "\t" +
                                        primers[row_pointer][column_pointer + 3] + "\t" +
                                        primers[row_pointer][column_pointer + 4])
                                primer_end_dict[row_pointer] = primer_end_set
                                primer_middle_dict[row_pointer] = primer_middle_set
                                primer_dict[row_pointer] = primer_set
                                primer_end_set = primer_end_set.union(current_inner_end)
                                primer_middle_set = primer_middle_set.union(current_inner_middle)
                                primer_set = primer_set.union(set(dege_trans(primers[row_pointer][column_pointer]) +
                                                                  dege_trans(primers[row_pointer][column_pointer + 1])))
                                jdict[row_pointer] = column_pointer
                                row_pointer += 1
                                break

    return clique


def greedy_maximal_primers(primers, row_num, clique):  # primers_file: primers file；row_num：row number of primers,equal to the cluster number; clique: clique
    global primer_end_set, maximal_check, primer_middle_set, primer_set
    primer_end_set = set()
    primer_middle_set = set()
    primer_dict = {}
    jdict = {}
    primer_end_dict = {}  
    primer_middle_dict = {} 
    row_pointer = 0
    while row_pointer < row_num:
        # i is Row pointer；j is Column pointer
        column_pointer = 1
        # print(i,len(primers[i]))
        if len(primers[row_pointer]) <= 1:
            print("virus {} missing!".format(primers[row_pointer][0]))
            row_pointer += 1
        else:
            while column_pointer <= len(primers[row_pointer]) - step:
                if inner_dimer_check(primers[row_pointer][column_pointer],
                                     primers[row_pointer][column_pointer + 1]) or \
                        outer_dimer_check(primers[row_pointer][column_pointer],
                                          primers[row_pointer][column_pointer + 1]):
                    column_pointer += step
                    if column_pointer > len(primers[row_pointer]) - step:
                        clique.append(primers[row_pointer][0] + "NA")
                        print("virus {} missing!".format(primers[row_pointer][0]))
                        row_pointer += 1
                        break
                else:
                    if column_pointer <= len(primers[row_pointer]) - step:
                        print(row_pointer, primers[row_pointer][0], column_pointer)  # ,global_primer_set)
                        clique.append(primers[row_pointer][0] + "\t" + str(column_pointer) + "\t" +
                                      primers[row_pointer][column_pointer] + "\t" +
                                      primers[row_pointer][column_pointer + 1] + "\t" +
                                      primers[row_pointer][column_pointer + 2] + "\t" +
                                      primers[row_pointer][column_pointer + 3] + "\t" +
                                      primers[row_pointer][column_pointer + 4])
                        primer_end_dict[row_pointer] = primer_end_set
                        primer_middle_dict[row_pointer] = primer_middle_set
                        primer_dict[row_pointer] = primer_set
                        primer_end_set = primer_end_set.union(current_inner_end)
                        primer_middle_set = primer_middle_set.union(current_inner_middle)
                        primer_set = primer_set.union(set(dege_trans(primers[row_pointer][column_pointer]) +
                                                          dege_trans(primers[row_pointer][column_pointer + 1])))
                        jdict[row_pointer] = column_pointer
                        row_pointer += 1
                        break
    return clique


def remove_None_and_sortby_len(list_list):
    list_list = [list(filter(None, i.strip().split("\t"))) for i in list_list]
    list_list = list(sorted(list_list, key=len))
    return list_list


def filter_len(l_element): # useless
    if len(l_element) > 1:
        return l_element


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    method = options.maxi
    kmer = options.kmer
    step = options.step
    distance = options.dist
    
    primers_file = open(options.input, "r")
    primers = primers_file.readlines()
    # print(primers[1].split("\t"))
    row_num = len(primers)
    # print(primers[0])
    # print(row_num)
    primers = remove_None_and_sortby_len(primers)
    # print(primers[0], primers[1])
    if re.search("/",options.input):
        sort_dir = options.input.split("/")
        sort = '/'.join(sort_dir[:-1]) + "/sort." + sort_dir[-1]
    else:
        sort = "sort." + options.input
    with open(sort, "w") as f:
        for i in primers:
            str_i = '\t'.join(i)
            f.write(str_i + "\n")

    if method == "T":
        maximal_primer_list = []
        maximal_out = options.out
        greedy_maximal_primers(primers, row_num, maximal_primer_list)
        with open(maximal_out, "w") as f:
            f.write("#primer\tprimer_rank\tPCR_product\t \
                primer_match_by_blast [2 mismatch]\t \
                primer_position_by_degePrimer_results\n" +
                '\n'.join(maximal_primer_list))
    else:
        max_primer_list = []
        if re.search("/",options.out):
            maximum_dir = options.out.split("/")
            maximum_out = '/'.join(maximum_dir[:-1]) + \
		"/maximum." + maximum_dir[-1]
        else:
            maximum_out = "maximum." + options.out
        greedy_primers(primers, row_num, max_primer_list)
        if maximal_check:
            pass
        else:
            with open(maximum_out, "w") as out:
                out.write('\n'.join(max_primer_list))

