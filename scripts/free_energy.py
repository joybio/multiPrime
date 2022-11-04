#!/bin/python

__date__ = "2022-10-8"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

import math
import sys

"""
The MIT License (MIT)

Copyright (c) 2022 Junbo Yang <yang_junbo_hi@126.com> <1806389316@pku.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from itertools import product
import pandas as pd
from optparse import OptionParser


def argsParse():
    parser = OptionParser('Usage: %prog -i input] -o [output] \n'
                          'Options: -g [gini] -p [position]', version="%prog 0.0.3")
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file.')
    parser.add_option('-g', '--gini',
                      dest='gini',
                      default="H_bonds",
                      type="str",
                      help='method used to calculate delta G. [unidfied] or [H_bonds]. Default: H_bonds.\n'
                           'NN model: nearest neighbor model')
    parser.add_option('-p', '--position',
                      dest='position',
                      default="0",
                      type="int",
                      help='which column is sequence. Default: [0]')
    parser.add_option('-f', '--format',
                      dest='format',
                      help='Format of primer file: xls [dataframe], fa [fasta] or seq; default: xls. \n'
                           'fa: fasta format. \n seq: sequence format, e.g. ATGCTGATGCATCGT.')

    parser.add_option('-o', '--out',
                      dest='out',
                      default="delta_G.txt",
                      type='str',
                      help='output file.')

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
    return parser.parse_args()


# TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")
TRANS = str.maketrans("ATGC", "TACG")

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

# Allawi, H. T. & SantaLucia, J., Jr. (1997) Biochemistry 36, 10581–10594
# adjust = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}
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
symmetry_correction = 0.4
#############################################################################
base2bit = {"A": 0, "C": 1, "G": 2, "T": 3}

# In duplex DNA there are 10 such unique doublets.
# These are 5 - 3 AA=TT; AG=CT; AC=GT; GA=TC; GG=CC; TG=CA, CG, GC, AT, and TA
freedom_of_degree_37_table_unified = pd.DataFrame({"A": [-1.00, -1.44, -1.28, -0.88],
                                                   "C": [-1.45, -1.84, -2.17, -1.28],
                                                   "G": [-1.30, -2.24, -1.84, -1.44],
                                                   "T": [-0.58, -1.30, -1.45, -1.00]},
                                                  index=["A", "C", "G", "T"])


def degenerate_seq(sequence):
    seq = []
    cs = ""
    for s in sequence:
        if s not in degenerate_base:
            cs += s
        else:
            seq.append([cs + i for i in degenerate_base[s]])
            cs = ""
    if cs:
        seq.append([cs])
    return ("".join(i) for i in product(*seq))


def delta_G(sequence, gini):
    Delta_G_list = []
    Na = 50
    Delta_G = 0
    if gini == "unified":
        for seq in degenerate_seq(sequence):
            i = 0
            while i < len(seq) - 1:
                Delta_G += freedom_of_degree_37_table_unified.loc[seq[i + 1], seq[i]]
                i += 1
            term5 = seq[-2:]
            if term5 == "TA":
                Delta_G += adjust_initiation[seq[0]] + adjust_initiation[seq[-1]] + adjust_terminal_TA
            else:
                Delta_G += adjust_initiation[seq[0]] + adjust_initiation[seq[-1]]
            # adjust by concentration of Na+
            Delta_G -= (0.175 * math.log(Na / 1000, math.e) + 0.20) * len(seq)
            if seq == seq[::-1]:
                Delta_G += symmetry_correction
            Delta_G_list.append(Delta_G)

    elif gini == "H_bonds":
        for seq in degenerate_seq(sequence):
            for n in range(len(seq) - 1):
                i, j = base2bit[seq[n + 1]], base2bit[seq[n]]
                Delta_G += freedom_of_H_37_table[i][j] * H_bonds_number[i][j] + penalty_of_H_37_table[i][j]
            term5 = sequence[-2:]
            if term5 == "TA":
                Delta_G += adjust_initiation[seq[0]] + adjust_initiation[seq[-1]] + adjust_terminal_TA
            else:
                Delta_G += adjust_initiation[seq[0]] + adjust_initiation[seq[-1]]
            Delta_G -= (0.175 * math.log(Na / 1000, math.e) + 0.20) * len(seq)
            if seq == seq[::-1]:
                Delta_G += symmetry_correction
            Delta_G_list.append(Delta_G)

    return round(max(Delta_G_list), 2)


if __name__ == "__main__":
    (options, args) = argsParse()
    if options.gini != "unified" and options.gini != "H_bonds":
        print("Method is wrong！Method must be unified or H_bonds.")
        sys.exit(1)
    else:
        if options.format == 'seq':
            sequence = options.input
            Delta_G = delta_G(sequence, options.gini)
            with open(options.out, "w") as output:
                output.write(sequence + "\t" + str(Delta_G) + "\n")
            output.close()
            print(Delta_G)
        else:
            f = open(options.input, "r")
            if options.format == 'xls':
                with open(options.out, "w") as output:
                    for i in f:
                        line = i.strip()
                        i = i.strip().split("\t")
                        sequence = i[options.position]
                        Delta_G = delta_G(sequence, options.gini)
                        output.write(line + "\t" + str(Delta_G) + "\n")
                output.close()
            elif options.format == 'fa':
                with open(options.out, "w") as output:
                    for i in f:
                        if i.startswith(">"):
                            output.write(i.strip() + "\t")
                        else:
                            sequence = i.strip()
                            Delta_G = delta_G(sequence, options.gini)
                            output.write(sequence + "\t" + str(Delta_G) + "\n")
                output.close()
            else:
                print("Check your input!!!")
                sys.exit(1)
            f.close()

