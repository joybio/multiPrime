#!/bin/python
"""
Get the final primerset.
"""
__date__ = "2022-7-6"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import os
import sys
from optparse import OptionParser
import re
import math
from operator import mul
from functools import reduce
import pandas as pd
import threading
from Bio.Seq import Seq


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
                      Default: TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT. \n \
                      If your sequence already have adapters, set -a ",".')

    parser.add_option('-f', '--form',
                      dest='form',
                      default="fa",
                      type="str",
                      help='File format: fa or txt. The txt format should like: primer_F_ID[\t]primer_F['
                           '\t]primer_R_ID[\t]primer_R[\n].')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='out file')
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Input file must be specified !!!")
        sys.exit(1)
    elif options.form is None:
        parser.print_help()
        print("file format must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("No output file provided !!!")
        sys.exit(1)
    return parser.parse_args()


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


#####################################################################################
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
                        break
        expand_score = reduce(mul, [score_trans(x) for x in expand_seq])
    return expand_seq


def Penalty_points(length, GC, d1, d2):
    return math.log((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)), 10)


#####################################################################################
def current_end(primer):
    primer_extend = adaptor[0] + primer
    # primer_len = len(primer)
    end_seq = set()
    for a in range(10):
        end_seq_primary = set(dege_trans(primer_extend[-a - 5:]))
        end_seq = end_seq.union(end_seq_primary)
    return end_seq


#####################################################################################
def get_primer_set(primers):
    primer_set = set()
    global primer_dict
    with open(primers, "r") as f:
        for i in f:
            i = i.strip().split("\t")
            primer_set.add(i[1])
            primer_dict[i[1]] = i[0]
            if len(i) == 2:
                pass
            elif len(i) == 4:
                primer_set.add(i[3])
                primer_dict[i[3]] = i[2]
    return primer_set


#####################################################################################
def dimer_check(primer, primer_set):
    current_end_set = current_end(primer)
    global matrix, end_mem
    for end in current_end_set:
        current_matrix = pd.DataFrame(
            columns=["Primer_ID", "Primer seq", "Primer end", "Delta G", "Primer end length", "End (distance 1)",
                     "End (GC)", "Dimer-primer_ID", "Dimer-primer seq", "End (distance 2)", "Loss"])
        if end in end_mem:
            if not end_mem[end].empty:
                end_mem[end]["Primer_ID"] = primer_dict[primer]
                end_mem[end]["Primer seq"] = primer
                matrix = pd.concat([matrix, end_mem[end]], axis=0, ignore_index=True)
        else:
            for ps in primer_set:
                expand_p = dege_trans(ps)
                for p in expand_p:
                    if re.search(str(Seq(end).reverse_complement()), p):
                        end_length = len(end)
                        end_GC = end.count("G") + end.count("C")
                        end_d1 = 0
                        end_d2 = len(p) - len(end) - p.index(str(Seq(end).reverse_complement()))
                        Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                        delta_G = delta_G_check(end)
                        if Loss > 3 or delta_G < -5:
                            matrix_local = pd.DataFrame({
                                "Primer_ID": primer_dict[primer],
                                "Primer seq": primer,
                                "Primer end": end,
                                "Delta G": delta_G,
                                "Primer end length": end_length,
                                "End (distance 1)": end_d1,
                                "End (GC)": end_GC,
                                "Dimer-primer_ID": primer_dict[ps],
                                "Dimer-primer seq": ps,
                                "End (distance 2)": end_d2,
                                "Loss": Loss
                            }, index=[0])
                            current_matrix = pd.concat([current_matrix, matrix_local], axis=0, ignore_index=True)
                            matrix = pd.concat([matrix, matrix_local], axis=0, ignore_index=True)
                    else:
                        pass
            end_mem[end] = current_matrix


#####################################################################################
def file_format(fasta, txt):
    out = open(txt, "w")
    n = 1
    with open(fasta, "r") as f:
        for i in f:
            if n%4 == 0:
                out.write(i)
            else:
                out.write(i.strip() + "\t")
            n += 1
    out.close()


#####################################################################################
freedom_of_H_37_table = pd.DataFrame({"A": [-0.70, -0.67, -0.69, -0.61],
                                      "C": [-0.81, -0.72, -0.87, -0.69],
                                      "G": [-0.65, -0.80, -0.72, -0.67],
                                      "T": [-0.65, -0.65, -0.81, -0.70]},
                                     index=["A", "C", "G", "T"])

penalty_of_H_37_table = pd.DataFrame({"A": [0.4, 0.23, 0.41, 0.33],
                                      "C": [0.575, 0.32, 0.45, 0.41],
                                      "G": [0.33, 0.17, 0.32, 0.23],
                                      "T": [0.73, 0.33, 0.575, 0.4]},
                                     index=["A", "C", "G", "T"])

H_bonds_number = pd.DataFrame({"A": [2, 2.5, 2.5, 2],
                               "C": [2.5, 3, 3, 2.5],
                               "G": [2.5, 3, 3, 2.5],
                               "T": [2, 2.5, 2.5, 2]},
                              index=["A", "C", "G", "T"])

adjust = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}


#####################################################################################
def delta_G_check(sequence):
    # 1   0%    25%    50%    75%   100%
    # 2 -13.61  -9.94  -9.06  -8.24  -5.90
    i = 0
    expand_seq = dege_trans(sequence[-9:])
    Delta_G_list = []
    for seq in expand_seq:
        Delta_G = 0
        while i < len(seq) - 1:
            Delta_G += (freedom_of_H_37_table.loc[seq[i + 1], seq[i]] * H_bonds_number.loc[
                seq[i + 1], seq[i]]) + \
                       penalty_of_H_37_table.loc[seq[i + 1], seq[i]]
            i += 1
        start = seq[0]
        stop = seq[-1]
        Delta_G += adjust[start] + adjust[stop]
        Delta_G_list.append(Delta_G)
    return round(max(Delta_G_list), 2)


#####################################################################################
if __name__ == "__main__":
    (options, args) = argsParse()
    adaptor = options.adaptor.split(",")
    adaptor_len = len(adaptor[0])
    if options.form == "fa":
        primers_file = options.input.rstrip(".fa") + ".txt"
        file_format(options.input, primers_file)
    elif options.form == "txt":
        primers_file = options.input
    else:
        print("Please check you input file format. Only fa or txt is accepted.")
        sys.exit(1)
    primer_dict = {}
    primer_set = get_primer_set(primers_file)
    end_mem = {}
    output = open(options.out, "w")
    matrix = pd.DataFrame(
        columns=["Primer_ID", "Primer seq", "Primer end", "Delta G", "Primer end length", "End (distance 1)",
                 "End (GC)",
                 "Dimer-primer_ID", "Dimer-primer seq",
                 "End (distance 2)", "Loss"])
    primers = open(primers_file, "r")
    for i in primers:
        i = i.strip().split("\t")
        if len(i) == 2:
            t = [threading.Thread(target=dimer_check, args=(i[1], primer_set))]
        elif len(i) == 4:
            t = [threading.Thread(target=dimer_check, args=(i[1], primer_set)),
                 threading.Thread(target=dimer_check, args=(i[3], primer_set))]
        else:
            break
        for t1 in t:
            t1.start()
            t1.join()
    matrix.to_csv(output, index=False, sep="\t")
    output.close()
    tmp = options.out + ".dimer_number"
    with open(tmp, "w") as t:
        primer_list = ["Primer_ID", "Dimer-primer_ID"]
        pieces = []
        for col in primer_list:
            tmp_series = matrix[col].value_counts()
            tmp_series.name = col
            pieces.append(tmp_series)
        df_value_counts = pd.concat(pieces, axis=1)
        df_value_counts.fillna(0, inplace=True)
        df_value_counts["rowSum"] = df_value_counts.apply(lambda x: x.sum(), axis=1)
        df_value_counts.sort_values(["rowSum"], ascending=False)
        df_value_counts.to_csv(t, index=True, sep="\t")
