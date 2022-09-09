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
                  help='File format: fa or txt. The txt format should like: primer_F_ID[\t]primer_F[\t]primer_R_ID[\t]primer_R[\n].')

parser.add_option('-o', '--out',
                  dest='out',
                  help='Prefix of out file: candidate primers')
(options, args) = parser.parse_args()

import re
import math
from operator import mul
from functools import reduce
import pandas as pd
import numpy as np
import threading

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


def current_end(primer):
    primer_extend = adaptor[0] + primer
    primer_len = len(primer)
    end_seq = set()
    for a in range(primer_len - 5):
        end_seq_primary = set(dege_trans(primer_extend[-a - 5:]))
        end_seq = end_seq.union(end_seq_primary)
        a += 1
    # print(end_seq)
    return end_seq


def get_primer_set(primers):
    primer_set = set()
    global primer_dict
    with open(primers,"r") as f:
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


def dimer_check(primer, primer_set):
    current_end_set = current_end(primer)
    global matrix
    # print(current_end_set)
    for p in primer_set:
        for end in current_end_set:
            if re.search(str(Seq(end).reverse_complement()), p):
                end_length = len(end)
                end_GC = end.count("G") + end.count("C")
                end_d1 = 0
                # print(end,p)
                end_d2 = len(p) - len(end) - p.index(str(Seq(end).reverse_complement()))
                Loss = Penalty_points(end_length, end_GC, end_d1, end_d2)
                if Loss > 3:
                    matrix_local = pd.DataFrame({
                        "Primer_ID": primer_dict[primer],
                        "Primer seq": primer,
                        "Primer end": end,
                        "Primer end length": end_length,
                        "End (distance 1)": end_d1,
                        "End (GC)": end_GC,
                        "Dimer-primer_ID": primer_dict[p],
                        "Dimer-primer seq": p,
                        "End (distance 2)": end_d2,
                        "Loss": Loss
                    }, index=[0])
                    matrix = pd.concat([matrix, matrix_local], axis=0, ignore_index=True)
                    # print(matrix)
            else:
                pass


def file_format(fasta, txt):
    n = 1
    out = open(txt, "w")
    with open(fasta, "r") as f:
        for i in f:
            if n % 4 == 0:
                out.write(i)
            else:
                out.write(i.strip() + "\t")
            n += 1
    out.close()


def remove_None_and_sortby_len(list_list):
    list_list = [list(filter(None, i.strip().split("\t"))) for i in list_list]
    list_list = list(sorted(list_list, key=len))
    return list_list


if __name__ == "__main__":
    adaptor = options.adaptor.split(",")
    adaptor_len = len(adaptor[0])
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

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
    output = open(options.out, "w")
    matrix = pd.DataFrame(
        columns=["Primer_ID", "Primer seq", "Primer end", "Primer end length", "End (distance 1)", "End (GC)", "Dimer-primer_ID", "Dimer-primer seq",
                 "End (distance 2)", "Loss"])
    primers = open(primers_file, "r")
    for i in primers:
        i = i.strip().split("\t")
        t = []
        if len(i) == 2:
            t.append(threading.Thread(target=dimer_check, args=(i[1], primer_set)))
        elif len(i) == 4:
            t.append(threading.Thread(target=dimer_check, args=(i[1], primer_set)))
            t.append(threading.Thread(target=dimer_check, args=(i[3], primer_set)))
        for t1 in t:
            t1.start()
        for t1 in t:
            t1.join()
    matrix.to_csv(output, index=False, sep="\t")
    tmp = options.out + ".dimer_number"
    with open(tmp,"w") as t:
        primer_count = pd.DataFrame(matrix.value_counts("Primer_ID"))
        primer_count.columns = ["number"]
        primer_count.sort_values(by="number",ascending=False)
        primer_count.to_csv(t,index=True,sep="\t")
        
    output.close()
    
    
