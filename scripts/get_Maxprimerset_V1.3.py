#!/bin/python
"""
Get the final primerset.
"""
__date__ = "2022-7-6"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import sys
from itertools import product
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
                      help='Adaptor sequence, which is used for NGS next (if not for NGS next, forget this '
                           'parameter). If you are not sure which adaptor will be used, You can use the '
                           'example sequence for dimer detection, because adaptor sequence will not form dimer '
                           'with primers generally. For example: TCTTTCCCTACACGACGCTCTTCCGATCT,'
                           'TCTTTCCCTACACGACGCTCTTCCGATCT. Default: TCTTTCCCTACACGACGCTCTTCCGATCT,'
                           'TCTTTCCCTACACGACGCTCTTCCGATCT. ')

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


def current_end(primer, num=5):
    end_seq = []
    for a in range(num, len(primer)):
        tmp_end_seq = dege_trans(primer[-a:])
        end_seq.extend(tmp_end_seq)
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


def dimer_examination(primer_F, primer_R, primer_set):
    def check_primer(primer, current_end_list):
        for end in current_end_list:
            idx = primer.find(RC(end))
            if idx >= 0:
                end_length = len(end)
                end_GC = end.count("G") + end.count("C")
                end_d1 = 0
                end_d2 = len(primer) - len(end) - idx
                Loss = Penalty_points(
                    end_length, end_GC, end_d1, end_d2)
                delta_G = deltaG(end)
                if Loss >= 3.96 or (delta_G < -5 and (end_d1 == end_d2)):
                    return True
        return False

    current_primer = set(degenerate_seq(primer_F) + degenerate_seq(primer_R))
    total_p_set = current_primer.union(primer_set)

    # 判断是否存在二聚体
    for p in current_primer:
        if check_primer(p, total_p_set):
            return True
    return False


def greedy_primers(primers, row_num, output):
    primer_end_set, primer_set = set(), set()
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    blank_row = 0
    column_pointer = 1
    # 构造空的dataframe，用于存储结果
    columns = ["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product (Length:Tm:Coverage)",
               "Coverage number with error in top N", "Primer position (representative sequence)"]
    clique = pd.DataFrame(columns=columns)

    def add_row_to_clique(row):
        nonlocal primer_end_set, primer_set, primer_end_dict, primer_dict, jdict, clique
        # 添加行
        clique_local = pd.DataFrame(row, index=[0], columns=columns)
        clique = pd.concat([clique, clique_local], axis=0, ignore_index=True)

        primer_end_dict[row_pointer] = primer_end_set
        primer_dict[row_pointer] = primer_set

        # 更新primer_end_set和primer_set
        current_end_set = current_end(row["Primer_F"]).union(current_end(row["Primer_R"]))
        primer_end_set = primer_end_set.union(current_end_set)
        primer_set = primer_set.union(set(dege_trans(row["Primer_F"]) + dege_trans(row["Primer_R"])))

        jdict[row_pointer] = column_pointer

    def backtrack_to_previous_row():
        nonlocal row_pointer, column_pointer, primer_end_set, primer_set, clique
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

    for row_pointer in range(row_num):
        if len(primers[row_pointer]) <= 1:
            blank_row += 1
        else:
            while column_pointer <= len(primers[row_pointer]) - step:
                if dimer_examination(primers[row_pointer][column_pointer],
                                     primers[row_pointer][column_pointer + 1], primer_set):
                    column_pointer += step
                    backtrack_to_previous_row()
                else:
                    add_row_to_clique({
                        "#Primer": primers[row_pointer][0],
                        "Primer_rank": str(column_pointer),
                        "Primer_F": primers[row_pointer][column_pointer],
                        "Primer_R": primers[row_pointer][column_pointer + 1],
                        "PCR_product (Length:Tm:Coverage)": primers[row_pointer][column_pointer + 2],
                        "Coverage number with error in top N": primers[row_pointer][column_pointer + 3],
                        "Primer position (representative sequence)": primers[row_pointer][column_pointer + 4]
                    })
                    column_pointer = 1
                    break

    # 结果输出到文件
    clique.to_csv(output, index=False, sep="\t")


def nan_removing(pre_list):
    while np.nan in pre_list:
        pre_list.remove(np.nan)
    return pre_list


def greedy_maximal_primers(primers, row_num, output, next_candidate):
    primer_end_set, primer_set = set(), set()
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    row_pointer = 0
    column_pointer = 1
    # 构造空的dataframe，用于存储结果
    columns = ["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product (Length:Tm:Coverage)",
               "Coverage number with error in top N", "Primer position (representative sequence)"]
    clique = pd.DataFrame(columns=columns)

    def add_row_to_clique(row):
        nonlocal primer_end_set, primer_set, primer_end_dict, primer_dict, jdict, clique
        clique_local = pd.DataFrame(row, index=[0], columns=columns)
        clique = pd.concat([clique, clique_local], axis=0, ignore_index=True)
        primer_end_dict[row_pointer] = primer_end_set
        primer_dict[row_pointer] = primer_set
        # 更新primer_end_set和primer_set
        current_end_set = current_end(row["Primer_F"]).union(current_end(row["Primer_R"]))
        primer_end_set = primer_end_set.union(current_end_set)
        primer_set = primer_set.union(set(dege_trans(row["Primer_F"]) + dege_trans(row["Primer_R"])))
        jdict[row_pointer] = column_pointer

    def process_normal_case():
        nonlocal row_pointer, column_pointer
        add_row_to_clique({
            "#Primer": primers[row_pointer][0],
            "Primer_rank": str(column_pointer),
            "Primer_F": primers[row_pointer][column_pointer],
            "Primer_R": primers[row_pointer][column_pointer + 1],
            "PCR_product (Length:Tm:Coverage)": primers[row_pointer][column_pointer + 2],
            "Coverage number with error in top N": primers[row_pointer][column_pointer + 3],
            "Primer position (representative sequence)": primers[row_pointer][column_pointer + 4]
        })
        row_pointer += 1
        column_pointer = 1

    def process_non_primer_case():
        nonlocal row_pointer, column_pointer
        next_virus = nan_removing(list(primers[row_pointer]))
        next_candidate.write("\t".join(next_virus) + "\n")
        row_pointer += 1
        column_pointer = 1

    while row_pointer < row_num:
        if len(primers[row_pointer]) <= 1:
            print("Non primers: virus {} missing!".format(primers[row_pointer][0]))
            process_non_primer_case()
        else:
            while column_pointer <= len(primers[row_pointer]) - step:
                if dimer_examination(primers[row_pointer][column_pointer], primers[row_pointer][column_pointer + 1],
                                     primer_set):
                    column_pointer += step
                    if column_pointer > len(primers[row_pointer]) - step:
                        clique_local = pd.DataFrame({"#Primer": primers[row_pointer][0]}, index=[0])
                        clique = pd.concat([clique, clique_local], axis=0, ignore_index=True)
                        print("virus {} missing!".format(primers[row_pointer][0]))
                        process_non_primer_case()
                        break
                else:
                    process_normal_case()
                    break

    # 结果输出到文件
    clique.to_csv(output, index=False, sep="\t")


if __name__ == "__main__":
    (options, args) = argsParse()
    adaptor = options.adaptor.split(",")
    adaptor_len = len(adaptor[0])
    if re.search("/", options.input):
        sort_dir = options.input.split("/")
        sort = '/'.join(sort_dir[:-1]) + "/sort." + sort_dir[-1]
    else:
        sort = "sort." + options.input
    with open(options.input, "r") as primers_file, open(sort, "w") as f:
        primers = list(sorted([list(filter(None, line.strip().split('\t'))) for line in primers_file], key=len))
        for i in primers:
            f.write('\t'.join(i) + "\n")
    method, step, row_num = options.method, options.step, len(primers)
    if method == "T":
        maximal_out = options.out
        next_candidate = options.out.rstrip(".xls") + ".next.xls"
        with open(next_candidate, "w") as next_candidate_txt:
            greedy_maximal_primers(primers, row_num, maximal_out, next_candidate_txt)
    else:
        maximum_out = options.out
        greedy_primers(primers, row_num, maximum_out)
