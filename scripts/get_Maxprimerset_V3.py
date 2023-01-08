#!/bin/python
"""
Get the final primerset.
construct primer_set by dimer_check, off_targets_prediction, deltaG_filter.
"""
__date__ = "2022-9-15"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import re
import os
import sys
from collections import defaultdict
from itertools import product

from Bio.Seq import Seq
from optparse import OptionParser
import re
import math
from operator import mul
from functools import reduce
import pandas as pd
import numpy as np
import threading
from pathlib import Path


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -o [output] \n \
                            Options: {-s [step] -m [T]}', version="%prog 0.0.4")
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: primers')

    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='reference file: bowtie index.')

    parser.add_option('-p', '--product',
                      dest='product',
                      default="150,2000",
                      help='Length of PCR product, default: [150,2000].')

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

    parser.add_option('-d', '--dist',
                      dest='dist',
                      default=12,
                      type="int",
                      help='{dist to 3end} primer length was used for mapping. Default: 12. '
                           'If -d 0, then the full primer will used for mapping.')

    parser.add_option('-t', '--tmp',
                      dest='tmp',
                      default="tmp",
                      type="str",
                      help='Temporary directory. Default: tmp.')

    parser.add_option('-m', '--method',
                      dest='method',
                      default="T",
                      type="str",
                      help='which method: maximal or maximum. Default: T. If -m [T] use maximal; else maximum')

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
    elif options.ref is None:
        parser.print_help()
        print("reference file (bowite index) must be specified !!!")
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

TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")

degenerate_base = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

score_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
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

TRANS_c = str.maketrans("ATCG", "TAGC")

def complement(seq):
    return seq.translate(TRANS_c)[::-1]
################################################################
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


################################################################
def Penalty_points(length, GC, d1, d2):
    return math.log((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)), 10)


def current_end(primer_F, primer_R):
    primer_F_extend = adaptor[0] + primer_F
    primer_R_extend = adaptor[1] + primer_R
    primer_F_len = len(primer_F)
    primer_R_len = len(primer_R)
    end_seq = set()
    for a in range(primer_F_len - 5):
        F_end_seq = set(dege_trans(primer_F_extend[-a - 5:]))
        end_seq = end_seq.union(F_end_seq)
        a += 1
    for b in range(primer_R_len - 5):
        R_end_seq = set(dege_trans(primer_R_extend[-b - 5:]))
        end_seq = end_seq.union(R_end_seq)
        b += 1
    return end_seq


################################################################
def dimer_check(primer_F, primer_R):
    check = "F"
    global primer_set
    current_end_set = current_end(primer_F, primer_R)
    for end in current_end_set:
        for p in primer_set:
            expand_p = dege_trans(p)
            for primer in expand_p:
                if re.search(str(Seq(end).reverse_complement()), primer):
                    end_length = len(end)
                    end_GC = end.count("G") + end.count("C")
                    end_d1 = 0
                    # if re.search(str(Seq(end).reverse_complement()), primer):
                    end_d2 = len(primer) - len(end) - primer.index(str(Seq(end).reverse_complement()))
                    Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                    delta_G = delta_G_check(end)
                else:
                    Loss = 0
                    delta_G = 0
                if Loss > 3 or delta_G < -9.94:
                    check = "T"
                    break
            if check == "T":
                break
        if check == "T":
            break
    if check == "T":
        return True
    else:
        return False

def degenerate_seq(primer):
    seq = []
    cs = ""
    for s in primer:
        if s not in degenerate_base:
            cs += s
        else:
            seq.append([cs + i for i in degenerate_base[s]])
            cs = ""
    if cs:
        seq.append([cs])
    return ("".join(i) for i in product(*seq))

def symmetry(seq):
    if len(seq) % 2 == 1:
        return False
    else:
        F = seq[:int(len(seq) / 2)]
        R = complement(seq[int(len(seq) / 2):][::-1])
        if F == R:
            return True
        else:
            return False

def delta_G_check(sequence):
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


################################################################


###############################################################
def greedy_primers(primers, row_num, output):
    global primer_end_set, primer_set
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    row_pointer = 0
    blank_row = 0
    column_pointer = 1
    F_dict, R_dict = defaultdict(list), defaultdict(list)
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
                    ID = str(Path(primers[row_pointer][0])).rstrip(".candidate.primers.txt")
                    # WAGTAGCCGAGAGATGCC GTAAGTTCATTGCCCACY 450 85 7100:7550
                    position = primers[row_pointer][column_pointer + 4].split(":")
                    F_seq = primers[row_pointer][column_pointer]
                    R_seq = primers[row_pointer][column_pointer + 1]
                    current_forward_dict, current_reverse_dict = term_map(ID, position, F_seq, R_seq)
                    F_dict, R_dict, tmp_forward_dict, tmp_reverse_dict = \
                        update_dict(current_forward_dict, current_reverse_dict, F_dict, R_dict)
                    if off_targets_prediction(tmp_forward_dict, tmp_reverse_dict):
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
                        F_dict.update(tmp_forward_dict)
                        R_dict.update(tmp_reverse_dict)
                        jdict[row_pointer] = column_pointer
                        row_pointer += 1
                        column_pointer = 1
                        break
    with open(output, "w") as greedy_primers_out:
        clique.to_csv(greedy_primers_out, index=False, sep="\t")


###################################################################
def nan_removing(pre_list):
    while np.nan in pre_list:
        pre_list.remove(np.nan)
    return pre_list


###################################################################
def get_term(ID, seq, dist, out):
    global term_bp  # 标记已经比对过的序列信息，再出现该序列信息则直接略过。
    with open(out, "w") as o:
        expand_seq = dege_trans(seq[-dist:])
        for seq_rank in range(len(expand_seq)):
            tmp_ID = ID + "_" + str(seq_rank)
            sequence = expand_seq[seq_rank]
            if sequence in term_bp:
                pass
            else:
                o.write(tmp_ID + "\n" + expand_seq[seq_rank] + "\n")
                term_bp.add(seq_rank)


###################################################################
def bowtie2_map(fa, ref_index, out, for_out, rev_out):
    os.system("bowtie2 -p 10 -f -N 1 -L 8 -a -x {} -f -U {} -S {}".format(ref_index, fa, out))
    # os.system("bowtie -f -v 0 -n 1 -a -p {} --best --strata {} {} -S {}".format(self.nproc, ref_index, fa, out))
    os.system("samtools view -F 16 {} > {}".format(out, for_out))
    os.system("samtools view -f 16 {} > {}".format(out, rev_out))


###################################################################
def term_map(ID, position, F_seq, R_seq):
    F_ID = ID + str(position[0]) + "_F"
    R_ID = ID + str(position[1]) + "_R"
    term9_F = Path(temp).joinpath(Path(F_ID)).with_stem("fa")
    term9_R = Path(temp).joinpath(Path(R_ID)).with_stem("fa")
    F_sam = Path(temp).joinpath(Path(F_ID)).with_stem("sam")
    F_for_sam = Path(temp).joinpath(Path(F_ID)).with_stem("for.sam")
    F_rev_sam = Path(temp).joinpath(Path(F_ID)).with_stem("rev.sam")
    R_sam = Path(temp).joinpath(Path(R_ID)).with_stem("sam")
    R_for_sam = Path(temp).joinpath(Path(R_ID)).with_stem("for.sam")
    R_rev_sam = Path(temp).joinpath(Path(R_ID)).with_stem("rev.sam")
    # current_forward_json = Path(temp).joinpath(Path(ID+"_forward")).with_stem("json")
    # current_reverse_json = Path(temp).joinpath(Path(ID+"_reverse")).with_stem("json")
    current_forward_dict, current_reverse_dict = defaultdict(list), defaultdict(list)
    if Path(term9_F).exists() and Path(term9_R).exists():
        pass
    else:
        if Path(term9_F).exists() and not Path(term9_R).exists():
            t1 = threading.Thread(target=get_term, args=(R_ID, R_seq, options.dist, term9_R))
            t1.start()
            t1.join()
            t2 = threading.Thread(target=bowtie2_map, args=(term9_R, path_2_ref, R_sam, R_for_sam, R_rev_sam))
            t2.start()
            t2.join()
            t3 = [threading.Thread(target=build_dict, args=(R_for_sam, current_forward_dict)),
                  threading.Thread(target=build_dict, args=(R_rev_sam, current_reverse_dict))]
            for t in t3:
                t.start()
            for t in t3:
                t.join()
        elif Path(term9_R).exists() and not Path(term9_F).exists():
            t1 = threading.Thread(target=get_term, args=(F_ID, F_seq, options.dist, term9_F))
            t1.start()
            t1.join()
            t2 = threading.Thread(target=bowtie2_map, args=(term9_F, path_2_ref, F_sam, F_for_sam, F_rev_sam))
            t2.start()
            t2.join()
            t3 = [threading.Thread(target=build_dict, args=(F_for_sam, current_forward_dict)),
                  threading.Thread(target=build_dict, args=(F_rev_sam, current_reverse_dict))]
            for t in t3:
                t.start()
            for t in t3:
                t.join()
        else:
            t1 = [threading.Thread(target=get_term, args=(F_ID, F_seq, options.dist, term9_F)),
                  threading.Thread(target=get_term, args=(R_ID, R_seq, options.dist, term9_R))]
            for t in t1:
                t.start()
            for t in t1:
                t.join()
            t2 = [threading.Thread(target=bowtie2_map, args=(term9_F, path_2_ref, F_sam, F_for_sam, F_rev_sam)),
                  threading.Thread(target=bowtie2_map, args=(term9_F, path_2_ref, R_sam, R_for_sam, R_rev_sam))]
            for t in t2:
                t.start()
            for t in t2:
                t.join()
            for_sam = Path(temp).joinpath(Path(F_ID + R_ID)).with_stem("for.sam")
            rev_sam = Path(temp).joinpath(Path(F_ID + R_ID)).with_stem("rev.sam")
            os.system("cat {} {} > {}".format(F_for_sam, R_for_sam, for_sam))
            os.system("cat {} {} > {}".format(F_rev_sam, R_rev_sam, rev_sam))
            t3 = [threading.Thread(target=build_dict, args=(for_sam, current_forward_dict)),
                  threading.Thread(target=build_dict, args=(rev_sam, current_reverse_dict))]
            for t1 in t3:
                t1.start()
            for t1 in t3:
                t1.join()
    return current_forward_dict, current_reverse_dict


#######################################################################
def build_dict(Input, current_dict):
    with open(Input, "r") as f:
        for i in f:
            i = i.strip().split("\t")
            primer = i[0]
            gene = i[2]
            primer_seq = i[1]
            primer_match_start = int(i[3]) - 1
            # current_dict[gene].append([primer, primer_match_start])
            gen_info = (primer, primer_match_start)
            current_dict[gene].append(gen_info)


#######################################################################
def expand_dict(current_dict, raw_dict, tmp_dict):
    for key in current_dict.keys():
        if key in raw_dict.keys():
            tmp_dict[key] = raw_dict[key].extend(current_dict[key])
        else:
            tmp_dict[key] = current_dict[key]


#######################################################################
def update_dict(current_forward_dict, current_reverse_dict, F_dict, R_dict):
    tmp_forward_dict = {}
    tmp_reverse_dict = {}
    if not F_dict and not R_dict:
        tmp_forward_dict = current_forward_dict
        tmp_reverse_dict = current_reverse_dict
    else:
        if not current_forward_dict and not current_reverse_dict:
            pass
        else:
            if current_forward_dict and not current_reverse_dict:
                t1 = threading.Thread(target=expand_dict, args=(current_forward_dict, F_dict, tmp_forward_dict))
                t1.start()
                t1.join()
            elif not current_forward_dict and current_reverse_dict:
                t1 = threading.Thread(target=expand_dict, args=(current_reverse_dict, R_dict, tmp_reverse_dict))
                t1.start()
                t1.join()
            else:
                t1 = [threading.Thread(target=expand_dict, args=(current_forward_dict, F_dict, tmp_forward_dict)),
                      threading.Thread(target=expand_dict, args=(current_reverse_dict, R_dict, tmp_reverse_dict))]
                for t in t1:
                    t.start()
                for t in t1:
                    t.join()
    return F_dict, R_dict, tmp_forward_dict, tmp_reverse_dict


def off_targets_prediction(tmp_forward_dict, tmp_reverse_dict):
    global prediction
    tmp_prediction = pd.DataFrame(
        columns=["Primer_F", "Primer_R", "Product length", "Chrom (or Genes)", "Start", "Stop"])
    for gene in tmp_forward_dict.keys():
        for primer_F in tmp_forward_dict[gene]:
            position_start = primer_F[1]
            for primer_R in tmp_reverse_dict[gene]:
                position_stop = primer_R[1]
                distance = int(position_stop) - int(position_start) + 1
                if int(product_len[0]) < distance < int(product_len[1]):
                    off_target = pd.DataFrame({"Primer_F": primer_F[0],
                                               "Primer_R": primer_R[0],
                                               "Product length": distance,
                                               "Chrom (or Genes)": gene,
                                               "Start": int(primer_F[1]) - 1,
                                               "Stop": int(primer_R[1]) - 1
                                               }, index=[0])
                    # 动态规划， 对比新加入的引物和off-target最大的引物
                    tmp_prediction = pd.concat([tmp_prediction, off_target], axis=0, ignore_index=True)
    primer_list = ["Primer_F", "Primer_R"]
    pieces = []
    for col in primer_list:
        tmp_series = tmp_prediction[col].value_counts()
        tmp_series.name = col
        pieces.append(tmp_series)
    df_value_counts = pd.concat(pieces, axis=1)
    df_value_counts.fillna(0, inplace=True)
    df_value_counts["rowSum"] = df_value_counts.apply(lambda x: x.sum(), axis=1)
    df_value_counts.sort_values(["rowSum"], ascending=False)
    if df_value_counts[df_value_counts["rowSum"] > 3000].empty():  # number of offtargets > 3000
        return True
    else:
        return False


def greedy_maximal_primers(primers, row_num, output, next_candidate):
    global primer_end_set, primer_set
    jdict = {}
    primer_end_dict = {}
    primer_dict = {}
    F_dict, R_dict = defaultdict(list), defaultdict(list)
    row_pointer = 0
    blank_row = 0
    column_pointer = 1
    clique = pd.DataFrame(columns=["#Primer", "Primer_rank", "Primer_F", "Primer_R", "PCR_product",
                                   "Primer target by blast+ [2 mismatch]", "Primer position (representative sequence)"])
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
                    # Cluster_0_87.candidate.primers.txt
                    ID = str(Path(primers[row_pointer][0])).rstrip(".candidate.primers.txt")
                    # WAGTAGCCGAGAGATGCC GTAAGTTCATTGCCCACY 450 85 7100:7550
                    position = primers[row_pointer][column_pointer + 4].split(":")
                    F_seq = primers[row_pointer][column_pointer]
                    R_seq = primers[row_pointer][column_pointer + 1]
                    current_forward_dict, current_reverse_dict = term_map(ID, position, F_seq, R_seq)
                    F_dict, R_dict, tmp_forward_dict, tmp_reverse_dict = \
                        update_dict(current_forward_dict, current_reverse_dict, F_dict, R_dict)
                    if off_targets_prediction(tmp_forward_dict, tmp_reverse_dict):
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
                        F_dict.update(tmp_forward_dict)
                        R_dict.update(tmp_reverse_dict)
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
    print("NGS adaptors: {} - {}".format(adaptor[0], adaptor[1]))
    reference = options.ref
    print("reference index (bowtie index): {}".format(reference))
    adaptor_len = len(adaptor[0])
    product_len = options.product.split(",")
    print("off-targets length: {} - {}".format(product_len[0], product_len[1]))
    path_2_ref = options.ref
    primer_end_set = set()
    primer_set = set()
    method = options.method
    step = options.step
    primers_file = open(options.input, "r")
    primers = primers_file.readlines()
    row_num = len(primers)
    primers = remove_None_and_sortby_len(primers)
    temp = Path(options.out).parent.joinpath(Path(options.tmp))
    os.system("mkdir {}".format(temp))
    prediction = pd.DataFrame(columns=["Primer_F", "Primer_R", "Product length", "Chrom (or Genes)", "Start", "Stop"])
    term_bp = set()
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
        greedy_maximal_primers(primers, row_num, maximal_out, temp)
        next_candidate_txt.close()
    else:
        maximum_out = options.out
        greedy_primers(primers, row_num, maximum_out)
