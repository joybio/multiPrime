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
from itertools import product
from sys import argv
from Bio.Seq import Seq
from optparse import OptionParser
import re
import math
from operator import mul
from functools import reduce
import pandas as pd
import numpy as np

def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -o [output] \n \
                            Options: {-s [step] -m [T]}', version="%prog 0.0.4")
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: primers')

    parser.add_option('-a', '--adaptor',
                      dest='adaptor',
                      default="TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT",
                      type="str",
                      help='Adaptor sequence, which is used for NGS next (if not for NGS next, forget this parameter). \n \
                      If you are not sure which adaptor will be used, \n \
                      You can use the example sequence for dimer detection,  \n \
                      because adaptor sequence will not form dimer with primers generally. \n \
                      For example: TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT. \n \
                      Default: TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT. ')

    parser.add_option('-s', '--step',
                      dest='step',
                      default=5,
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
    # while count:
    #     try:
    #         import Bio  #
    #
    #         print('Dependent package Biopython is OK.\nDependent module Bio is OK.')
    #         break
    #     except:
    #         print('Dependent package Biopython is not found!!! \n Start installing ....')
    #         os.system('pip install biopython')
    #         count -= 1
    #         continue
    return parser.parse_args()


degenerate_pair = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

degenerate_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
                    "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}

# Martin Zacharias* "Base-Pairing and Base-Stacking Contributions to Double-Stranded DNA Formation"
# J. Phys. Chem. B 2020, 124, 46, 10345–10352

freedom_of_H_37_table = [[-0.7, -0.81, -0.65, -0.65],
                         [-0.67, -0.72, -0.8, -0.65],
                         [-0.69, -0.87, -0.72, -0.81],
                         [-0.61, -0.69, -0.67, -0.7]]

penalty_of_H_37_table = [[0.4, 0.575, 0.33, 0.73],
                         [0.23, 0.32, 0.17, 0.33],
                         [0.41, 0.45, 0.32, 0.575],
                         [0.33, 0.41, 0.23, 0.4]]

H_bonds_number = [[2, 2.5, 2.5, 2],
                  [2.5, 3, 3, 2.5],
                  [2.5, 3, 3, 2.5],
                  [2, 2.5, 2.5, 2]]
#############################################################################
# Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability
# adjust_initiation = {"A": 2.8, "T": 2.8, "C": 1.82, "G": 1.82}
# deltaG(total) = Σ(deltaG(i)) + deltaG(initiation with terminal GC) + deltaG(initiation with terminal AT) -
# (0.175 * ln(Na) + 0.2) * len(sequence)
# This work suggested that oligonucleotides with terminal 5-T-A-3 base pairs should have a penalty of 0.4 kcal/mol
# but that no penalty should be given for terminal 5-A-T-3 pairs
adjust_initiation = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}

adjust_terminal_TA = 0.4
# Symmetry correction applies only to self-complementary sequences.
# symmetry_correction = 0.4
symmetry_correction = 0.4
#############################################################################
base2bit = {"A": 0, "C": 1, "G": 2, "T": 3}

TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")

def RC(seq):
    return seq.translate(TRANS)[::-1]

def score_trans(sequence):
    return reduce(mul, [math.floor(degenerate_table[x]) for x in list(sequence)])

def symmetry(seq):
    if len(seq) % 2 == 1:
        return False
    else:
        F = seq[:int(len(seq) / 2)]
        R = RC(seq[int(len(seq) / 2):][::-1])
        if F == R:
            return True
        else:
            return False

def dege_trans(primer):
    seq = []
    cs = ""
    for s in primer:
        if s not in degenerate_pair:
            cs += s
        else:
            seq.append([cs + i for i in degenerate_pair[s]])
            cs = ""
    if cs:
        seq.append([cs])
    return ["".join(i) for i in product(*seq)]

def Penalty_points(length, GC, d1, d2):
    # return log10((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)))
    return math.log10((2 ** length * 2 ** GC) / ((2 ** d1 - 0.9) * (2 ** d2 - 0.9)))


def current_end(primer_F, primer_R, num=5, length=14):
    primer_F_extend = adaptor[0] + primer_F
    primer_R_extend = adaptor[1] + primer_R
    end_seq = []
    for a in range(num, (num + length)):
        F_end_seq = dege_trans(primer_F_extend[-a:])
        end_seq.extend(F_end_seq)
    for b in range(num, (num + length)):
        R_end_seq = set(dege_trans(primer_R_extend[-b:]))
        end_seq.extend(R_end_seq)
    return set(end_seq)


def degenerate_seq(primer):
    seq = []
    cs = ""
    for s in primer:
        if s not in degenerate_pair:
            cs += s
        else:
            seq.append([cs + i for i in degenerate_pair[s]])
            cs = ""
    if cs:
        seq.append([cs])
    return ["".join(i) for i in product(*seq)]

def deltaG(sequence):
    Delta_G_list = []
    Na = 50
    for seq in degenerate_seq(sequence):
        Delta_G = 0
        for n in range(len(seq) - 1):
            i, j = base2bit[seq[n + 1]], base2bit[seq[n]]
            Delta_G += freedom_of_H_37_table[i][j] * H_bonds_number[i][j] + penalty_of_H_37_table[i][j]
        term5 = sequence[-2:]
        if term5 == "TA":
            Delta_G += adjust_initiation[seq[0]] + adjust_initiation[seq[-1]] + adjust_terminal_TA
        else:
            Delta_G += adjust_initiation[seq[0]] + adjust_initiation[seq[-1]]
        # adjust by concentration of Na+
        Delta_G -= (0.175 * math.log(Na / 1000, math.e) + 0.20) * len(seq)
        if symmetry(seq):
            Delta_G += symmetry_correction
        Delta_G_list.append(Delta_G)
    return round(max(Delta_G_list), 2)

def dimer_check(primer_F, primer_R):
    check = "F"
    global primer_set
    current_end_set = current_end(primer_F, primer_R)
    print(current_end_set)
    current_end_list = sorted(list(current_end_set), key=lambda i: len(i), reverse=True)
    for end in current_end_list:
        for primer in primer_set:
            idx = primer.find(RC(end))
            if idx >= 0:
                end_length = len(end)
                end_GC = end.count("G") + end.count("C")
                end_d1 = 0
                end_d2 = len(primer) - len(end) - idx
                Loss = Penalty_points(
                    end_length, end_GC, end_d1, end_d2)
                delta_G = deltaG(end)
                if Loss >= 3 or (delta_G < -5 and (end_d1 == end_d2)):
                    check = "T"
                    break
        if check == "T":
            break
    if check == "T":
        return True
    else:
        set_end_seq = set()
        for primer in primer_set:
            for a in range(5, 19):
                tmp_end_seq = set(dege_trans(primer[-a:]))
                set_end_seq = set_end_seq.union(tmp_end_seq)
        for end in set_end_seq:
            for primer in [primer_F, primer_R]:
                idx = primer.find(RC(end))
                if idx >= 0:
                    end_length = len(end)
                    end_GC = end.count("G") + end.count("C")
                    end_d1 = 0
                    end_d2 = len(primer) - len(end) - idx
                    Loss = Penalty_points(
                        end_length, end_GC, end_d1, end_d2)
                    delta_G = deltaG(end)
                    if Loss >= 3 or (delta_G < -5 and (end_d1 == end_d2)):
                        check = "T"
                        break
            if check == "T":
                break
        if check == "T":
            return True
        else:
            return False


###############################################################
def greedy_primers(primers, row_num, output):
    global primer_end_set, primer_set
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    row_pointer = 0
    blank_row = 0
    column_pointer = 1
    clique = pd.DataFrame(columns=["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product (Length:Tm:Coverage)",
                                   "Coverage number with error in top N", "Primer position (representative sequence)"])
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
                        "PCR_product (Length:Tm:Coverage)": primers[row_pointer][column_pointer + 2],
                        "Coverage number with error in top N": primers[row_pointer][column_pointer + 3],
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

def nan_removing(pre_list):
    while np.nan in pre_list:
        pre_list.remove(np.nan)
    return pre_list

def greedy_maximal_primers(primers, row_num, output,next_candidate):
    global primer_end_set, primer_set
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    row_pointer = 0
    blank_row = 0
    column_pointer = 1
    clique = pd.DataFrame(columns=["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product (Length:Tm:Coverage)",
                                   "Coverage number with error in top N", "Primer position (representative sequence)"])
    while row_pointer < row_num:
        if len(primers[row_pointer]) <= 1:
            print("Non primers: virus {} missing!".format(primers[row_pointer][0]))
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
                        next_virus = nan_removing(list(primers[row_pointer]))
                        next_candidate.write("\t".join(next_virus) + "\n")
                        row_pointer += 1
                        column_pointer = 1
                        break
                    
                else:
                    clique_local = pd.DataFrame({
                        "#Primer": primers[row_pointer][0],
                        "Primer_rank": str(column_pointer),
                        "Primer_F": primers[row_pointer][column_pointer],
                        "Primer_R": primers[row_pointer][column_pointer + 1],
                        "PCR_product (Length:Tm:Coverage)": primers[row_pointer][column_pointer + 2],
                        "Coverage number with error in top N": primers[row_pointer][column_pointer + 3],
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
    (options, args) = argsParse()
    adaptor = options.adaptor.split(",")
    adaptor_len = len(adaptor[0])
    primer_end_set = set()
    primer_set = set()
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
        next_candidate = options.out.rstrip(".xls") + ".next.xls"
        next_candidate_txt = open(next_candidate, "w")
        greedy_maximal_primers(primers, row_num, maximal_out,next_candidate_txt)
        next_candidate_txt.close()
    else:
        maximum_out = options.out
        greedy_primers(primers, row_num, maximum_out)
