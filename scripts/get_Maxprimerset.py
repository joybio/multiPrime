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
                        Options: {-s [step] -m [T]}', version="%prog 0.0.4")
parser.add_option('-i', '--input',
                  dest='input',
                  help='Input file: primers')

parser.add_option('-a', '--adaptor',
                  dest='adaptor',
                  default=",",
                  type="str",
                  help='Adaptor sequence, which is used for NGS next (if not for NGS next, forget this parameter). \n \
                  If you are not sure which adaptor will be used, \n \
                  You can use the example sequence for dimer detection,  \n \
                  because adaptor sequence will not form dimer with primers generally. \n \
                  For example: TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT. \n \
                  Default: None. ')

parser.add_option('-s', '--step',
                  dest='step',
                  default=4,
                  type="int",
                  help='distance between primers; column number of primer1_F to primer2_F.')

parser.add_option('-m', '--method',
                  dest='method',
                  default="T",
                  type="str",
                  help='which method: maximal or maximum. If -m [T] use maximal; else maximum')

parser.add_option('-o', '--out',
                  dest='out',
                  help='Prefix of out file: candidate primers')
(options, args) = parser.parse_args()

import re
import math
from operator import mul
from functools import reduce
import pandas as pd

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


def Penalty_points(length, GC, d1, d2):
    return math.log((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)), 10)


adaptor = options.adaptor.split(",")
adaptor_len = len(adaptor[0])


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


def dimer_check(primer_F, primer_R):
    check = "F"
    global primer_set
    current_end_set = current_end(primer_F, primer_R)
    for end in current_end_set:
        for primer in primer_set:
            if re.search(str(Seq(end).reverse_complement()), primer):
                end_length = len(end)
                end_GC = end.count("G") + end.count("C")
                end_d1 = 0
                if re.search(str(Seq(end).reverse_complement()), primer):
                    end_d2 = min((len(primer) - len(end) - primer.index(str(Seq(end).reverse_complement()))),
                                 primer.index(str(Seq(end).reverse_complement())) + adaptor_len)
                    Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                else:
                    Loss = 0
                if Loss > 3:
                    check = "T"
                    break
        if check == "T":
            break
    if check == "T":
        return True
    else:
        return False


###############################################################
primer_end_set = set()
primer_set = set()

def greedy_primers(primers, row_num, output):
    global primer_end_set, primer_set
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    row_pointer = 0
    blank_row = 0
    column_pointer = 1
    clique = pd.DataFrame(columns=["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product",
                                   "Primer target by blast+ [2 mismatch]", "Primer position (representative sequence)"])
    while row_pointer < row_num:
        if len(primers[row_pointer]) <= 1:
            row_pointer += 1
            blank_row += 1
        else:
            while column_pointer <= len(primers[row_pointer]) - step:
                if dimer_check(primers[row_pointer][column_pointer],
                               primers[row_pointer][column_pointer + 1]):
                    column_pointer += step
                    while column_pointer > len(primers[row_pointer]) - step:
                        row_pointer -= 1
                        if row_pointer < blank_row:
                            print("Non maximum primer set. Try maximal primer set!")
                            sys.exit(1)
                        else:
                            column_pointer = jdict[row_pointer] + step
                            primer_end_set = primer_end_dict[row_pointer]
                            primer_set = primer_dict[row_pointer]
                            clique.drop([len(clique) - 1], inplace=True)
                else:
                    clique_local = pd.DataFrame({
                        "#Primer": primers[row_pointer][0],
                        "Primer_rank": str(column_pointer),
                        "Primer_F": primers[row_pointer][column_pointer],
                        "Primer_R": primers[row_pointer][column_pointer + 1],
                        "PCR_product": primers[row_pointer][column_pointer + 2],
                        "Primer target number by blast+ [2 mismatch]": primers[row_pointer][column_pointer + 3],
                        "Primer position (representative sequence)": primers[row_pointer][column_pointer + 4]
                    }, index=[0])
                    clique = pd.concat([clique, clique_local], axis=0, ignore_index=True)
                    primer_end_dict[row_pointer] = primer_end_set
                    primer_dict[row_pointer] = primer_set
                    current_end_set = current_end(primers[row_pointer][column_pointer],
                                                  primers[row_pointer][column_pointer + 1])
                    primer_end_set = primer_end_set.union(current_end_set)
                    primer_set = primer_set.union(set(dege_trans(primers[row_pointer][column_pointer]) +
                                                      dege_trans(primers[row_pointer][column_pointer + 1])))
                    jdict[row_pointer] = column_pointer
                    row_pointer += 1
                    column_pointer = 1
                    break
    with open(output, "w") as greedy_primers_out:
        clique.to_csv(greedy_primers_out, index=False, sep="\t")


def greedy_maximal_primers(primers, row_num, output):
    global primer_end_set, primer_set
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    row_pointer = 0
    blank_row = 0
    column_pointer = 1
    clique = pd.DataFrame(columns=["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product",
                                   "Primer target by blast+ [2 mismatch]", "Primer position (representative sequence)"])
    while row_pointer < row_num:
        if len(primers[row_pointer]) <= 1:
            print("virus {} missing!".format(primers[row_pointer][0]))
            row_pointer += 1
            blank_row += 1
        else:
            while column_pointer <= len(primers[row_pointer]) - step:
                if dimer_check(primers[row_pointer][column_pointer],
                               primers[row_pointer][column_pointer + 1]):
                    column_pointer += step
                    if column_pointer > len(primers[row_pointer]) - step:
                        clique_local = pd.DataFrame({"#Primer": primers[row_pointer][0]}, index=[0])
                        clique = pd.concat([clique, clique_local], axis=0, ignore_index=True)
                        print("virus {} missing!".format(primers[row_pointer][0]))
                        row_pointer += 1
                else:
                    clique_local = pd.DataFrame({
                        "#Primer": primers[row_pointer][0],
                        "Primer_rank": str(column_pointer),
                        "Primer_F": primers[row_pointer][column_pointer],
                        "Primer_R": primers[row_pointer][column_pointer + 1],
                        "PCR_product": primers[row_pointer][column_pointer + 2],
                        "Primer target number by blast+ [2 mismatch]": primers[row_pointer][column_pointer + 3],
                        "Primer position (representative sequence)": primers[row_pointer][column_pointer + 4]
                    }, index=[0])
                    clique = pd.concat([clique, clique_local], axis=0, ignore_index=True)
                    primer_end_dict[row_pointer] = primer_end_set
                    primer_dict[row_pointer] = primer_set
                    current_end_set = current_end(primers[row_pointer][column_pointer],
                                                  primers[row_pointer][column_pointer + 1])
                    primer_end_set = primer_end_set.union(current_end_set)
                    primer_set = primer_set.union(set(dege_trans(primers[row_pointer][column_pointer]) +
                                                      dege_trans(primers[row_pointer][column_pointer + 1])))
                    jdict[row_pointer] = column_pointer
                    row_pointer += 1
                    column_pointer = 1
                    break
    with open(output, "w") as greedy_maximal_primers_out:
        clique.to_csv(greedy_maximal_primers_out, index=False, sep="\t")


def remove_None_and_sortby_len(list_list):
    list_list = [list(filter(None, i.strip().split("\t"))) for i in list_list]
    list_list = list(sorted(list_list, key=len))
    return list_list


def filter_len(l_element):  # useless
    if len(l_element) > 1:
        return l_element


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    method = options.method
    step = options.step
    primers_file = open(options.input, "r")
    primers = primers_file.readlines()
    row_num = len(primers)
    primers = remove_None_and_sortby_len(primers)
    if re.search("/", options.input):
        sort_dir = options.input.split("/")
        sort = '/'.join(sort_dir[:-1]) + "/sort." + sort_dir[-1]
    else:
        sort = "sort." + options.input
    with open(sort, "w") as f:
        for i in primers:
            str_i = '\t'.join(i)
            f.write(str_i + "\n")
    if method == "T":
        maximal_out = options.out
        greedy_maximal_primers(primers, row_num, maximal_out)
    else:
        maximum_out = options.out
        greedy_primers(primers, row_num, maximum_out)
