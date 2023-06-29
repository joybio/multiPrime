#!/bin/python
# Bug fixed. 
# Sometimes, positons have no bases except "-", and columns in frequency array become [0,0,0,0],
# which is not proper for primer design, especially for Tm calculation,
# we repalce "-" in the start and stop region with the flanking nucloetides.
# we also set H and S to 0 in this situation to continue Tm calculation.


__date__ = "2023-4-3"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

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

import math
import time
from functools import reduce
from math import log10
from itertools import product
from multiprocessing import Manager
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor
import argparse
from operator import mul
from statistics import mean
import sys
from itertools import repeat
import json
import numpy as np
import pandas as pd


# Melting temperature between 55-80◦C reduces the occurrence of hairpins
# Runs of three or more Cs or Gs at the 3'-ends of primers may promote mispriming at G or C-rich sequences
# (because of stability of annealing), and should be avoided.
def parseArg():
    parser = argparse.ArgumentParser(description="For degenerate primer design")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input file: multi-alignment output (muscle or others).", metavar="<file>")
    parser.add_argument("-l", "--plen", type=int, default=18,
                        help='Length of primer. Default: 18.', metavar="<int>")
    parser.add_argument("-n", "--dnum", type=int, default=4,
                        help='Max number of degenerate. Default: 4.', metavar="<int>")
    parser.add_argument("-d", "--degeneracy", type=int, default=10,
                        help='Max degeneracy of primer. Default: 10.', metavar="<int>")
    parser.add_argument("-v", "--variation", type=int, default=1,
                        help='Max mismatch number of primer. Default: 1', metavar="<int>")
    parser.add_argument("-e", "--entropy", type=float, default=3.6,
                        help='Entropy is actually a measure of disorder. This parameter is used to judge whether the '
                             'window is conservation. Entropy of primer-length window. Default: 3.6.',
                        metavar="<float>")
    parser.add_argument("-g", "--gc", type=str, default="0.2,0.7",
                        help='Filter primers by GC content. Default [0.2,0.7].', metavar="<str>")
    parser.add_argument("-s", "--size", type=int, default=100,
                        help='Filter primers by mini PRODUCT size. Default 100.', metavar="<int>")
    parser.add_argument("-f", "--fraction", type=float, default=0.8,
                        help='Filter primers by match fraction. If you set -s lower than 0.8, make sure that'
                             '--entropy greater than 3.6, because disorder region (entropy > 3.6) will not be processed'
                             'in multiPrime. Even these regions can design coverage with error greater than your '
                             'threshold, it wont be processed. Default: 0.8.', metavar="<float>")
    parser.add_argument("-c", "--coordinate", type=str, default="1,2,-1",
                        help='Mismatch index is not allowed to locate in your specific positions.'
                             'otherwise, it wont be regard as the mis-coverage. With this param, '
                             'you can control the index of Y-distance (number=variation and position of mismatch)'
                             'when calculate coverage with error. coordinate>0: 5\'==>3\'; coordinate<0: 3\'==>5\'.'
                             'You can set this param to any value that you prefer. Default: 1,-1. '
                             '1:  I dont want mismatch at the 2nd position, start from 0.'
                             '-1: I dont want mismatch at the -1st position, start fro -1.', metavar="<str>")
    parser.add_argument("-p", "--proc", type=int, default=20,
                        help='Number of process to launch. Default: 20.', metavar="<int>")
    parser.add_argument("-a", "--away", type=int, default=4,
                        help='Filter hairpin structure, which means distance of the minimal paired bases. Default: 4. '
                             'Example:(number of X) AGCT[XXXX]AGCT. '
                             'Primers should not have complementary sequences (no consecutive 4 bp complementarities),'
                             'otherwise the primers themselves will fold into hairpin structure.', metavar="<int>")
    parser.add_argument("-o", "--out", type=str, required=True,
                        help='output file', metavar="<file>")
    return parser.parse_args()


degenerate_base = {"-": ["-"], "A": ["A"], "G": ["G"], "C": ["C"], "T": ["T"], "R": ["A", "G"], "Y": ["C", "T"],
                   "M": ["A", "C"], "K": ["G", "T"], "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"],
                   "B": ["G", "T", "C"], "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

score_table = {"-": 100, "#": 0.00, "A": 1, "G": 1.11, "C": 1.21, "T": 1.40, "R": 2.11, "Y": 2.61, "M": 2.21,
               "K": 2.51, "S": 2.32, "W": 2.40, "H": 3.61, "B": 3.72, "V": 3.32, "D": 3.51, "N": 4.72}

trans_score_table = {v: k for k, v in score_table.items()}

##############################################################################################
############################# Calculate free energy ##########################################
##############################################################################################
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
adjust_initiation = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}
adjust_terminal_TA = 0.4
# Symmetry correction applies only to self-complementary sequences.
# symmetry_correction = 0.4
symmetry_correction = 0.4

##############################################################################################
base2bit = {"A": 0, "C": 1, "G": 2, "T": 3, "#": 4}

##############################################################################################
# 37°C and 1 M NaCl
Htable2 = [[-7.9, -8.5, -8.2, -7.2, 0],
           [-8.4, -8, -9.8, -8.2, 0],
           [-7.8, -10.6, -8, -8.5, 0],
           [-7.2, -7.8, -8.4, -7.9, 0],
           [0, 0, 0, 0, 0]]
Stable2 = [[-22.2, -22.7, -22.2, -21.3, 0],
           [-22.4, -19.9, -24.4, -22.2, 0],
           [-21, -27.2, -19.9, -22.7, 0],
           [-20.4, -21, -22.4, -22.2, 0],
           [0, 0, 0, 0, 0]]
Gtable2 = [[-1, -1.45, -1.3, -0.58, 0],
           [-1.44, -1.84, -2.24, -1.3, 0],
           [-1.28, -2.17, -1.84, -1.45, 0],
           [-0.88, -1.28, -1.44, -1, 0],
           [0, 0, 0, 0, 0]]
H_adjust_initiation = {"A": 2.3, "T": 2.3, "C": 0.1, "G": 0.1}
S_adjust_initiation = {"A": 4.1, "T": 4.1, "C": -2.8, "G": -2.8}
G_adjust_initiation = {"A": 1.03, "T": 1.03, "C": 0.98, "G": 0.98}
H_symmetry_correction = 0
S_symmetry_correction = -1.4
G_symmetry_correction = 0.4
##############################################################################################
# ng/ul
primer_concentration = 100
Mo_concentration = 50
Di_concentration = 1.5
dNTP_concentration = 0.25
Kelvin = 273.15
# reference (Owczarzy et al.,2008)
crossover_point = 0.22

bases = np.array(["A", "C", "G", "T"])
di_bases = []
for i in bases:
    for j in bases:
        di_bases.append(i + j)


def Penalty_points(length, GC, d1, d2):
    return log10((2 ** length * 2 ** GC) / ((2 ** d1 - 0.9) * (2 ** d2 - 0.9)))


di_nucleotides = set()
for i in base2bit.keys():
    single = i * 4
    di_nucleotides.add(single)
    for j in base2bit.keys():
        if i != j:
            di = (i + j) * 4
            di_nucleotides.add(di)
        for k in base2bit.keys():
            if i != j != k:
                tri = (i + j + k) * 3
                di_nucleotides.add(tri)


def score_trans(sequence):
    return reduce(mul, [math.floor(score_table[x]) for x in list(sequence)])


def dege_number(sequence):
    return sum(math.floor(score_table[x]) > 1 for x in list(sequence))


TRANS = str.maketrans("ATGC", "TACG")


def RC(seq):
    return seq.translate(TRANS)[::-1]


##############################################################################################
############## m_distance which is used to calculate (n)-nt variation coverage ###############
# Caution: this function only works when degeneracy of seq2 < 2 (no degenerate in seq2).
##############################################################################################
def Y_distance(seq1, seq2):
    seq_diff = list(np.array([score_table[x] for x in list(seq1)]) - np.array([score_table[x] for x in list(seq2)]))
    m_dist = [idx for idx in range(len(seq_diff)) if round(seq_diff[idx], 2) not in score_table.values()]
    # print(seq_diff)
    return m_dist


##############################################################################################
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


def Calc_deltaH_deltaS(seq):
    Delta_H = 0
    Delta_S = 0
    for n in range(len(seq) - 1):
        i, j = base2bit[seq[n + 1]], base2bit[seq[n]]
        Delta_H += Htable2[i][j]
        Delta_S += Stable2[i][j]
    seq = seq.replace("#", '')
    Delta_H += H_adjust_initiation[seq[0]] + H_adjust_initiation[seq[-1]]
    Delta_S += S_adjust_initiation[seq[0]] + S_adjust_initiation[seq[-1]]
    if symmetry(seq):
        Delta_S += S_symmetry_correction
    return Delta_H * 1000, Delta_S


# salt_adjust = math.log(Tm_Na_adjust / 1000.0, math.e)
# def S_adjust(seq):
#     n = len(seq) - 1
#     # S_Na_adjust = 0.847 * n * salt_adjust
#     # Oligonucleotide Melting Temperatures under PCR Conditions: Nearest-Neighbor Corrections for
#     # Mg2+ , Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide Concentrations with
#     # Comparison to Alternative Empirical Formulas
#     S_Na_adjust = 0.368 * n * salt_adjust
#     # A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics
#     return S_Na_adjust
# where n is the total number of phosphates in the duplex divided by 2,
# This is equal to the oligonucleotide length minus 1.

def GC_fraction(seq):
    return round((list(seq).count("G") + list(seq).count("C")) / len(list(seq)), 3)


# different salt corrections for monovalent (Owczarzy et al.,2004) and divalent cations (Owczarzy et al.,2008)
def Calc_Tm_v2(seq):
    delta_H, delta_S = Calc_deltaH_deltaS(seq)
    # Note that the concentrations in the following Eq is mmol/L, In all other equations,concentration are mol/L
    # Monovalent cations are typically present as K+ and Tris+ in PCR buffer,
    # K+ is similar to Na+ in regard to duplex stabilization
    # if Di_concentration > dNTP_concentration:
    #     Tm_Na_adjust = Mo_concentration + 120 * math.sqrt(Di_concentration - dNTP_concentration)
    # else:
    #     Tm_Na_adjust = Mo_concentration
    Tm_Na_adjust = Mo_concentration

    if dNTP_concentration >= Di_concentration:
        free_divalent = 0.00000000001
    else:
        free_divalent = (Di_concentration - dNTP_concentration) / 1000.0
    R_div_monov_ratio = (math.sqrt(free_divalent)) / (Mo_concentration / 1000)

    if R_div_monov_ratio < crossover_point:
        # use only monovalent salt correction, [equation 22] (Owczarzy et al., 2004)
        correction = (((4.29 * GC_fraction(seq)) - 3.95) * pow(10, -5) * math.log(Tm_Na_adjust / 1000.0, math.e)) \
                     + (9.40 * pow(10, -6) * (pow(math.log(Tm_Na_adjust / 1000.0, math.e), 2)))
    else:
        # magnesium effects are dominant, [equation 16] (Owczarzy et al., 2008) is used
        # Table 2
        a = 3.92 * pow(10, -5)
        b = - 9.11 * pow(10, -6)
        c = 6.26 * pow(10, -5)
        d = 1.42 * pow(10, -5)
        e = - 4.82 * pow(10, -4)
        f = 5.25 * pow(10, -4)
        g = 8.31 * pow(10, -5)
        if R_div_monov_ratio < 6.0:
            a = 3.92 * pow(10, -5) * (
                    0.843 - (0.352 * math.sqrt(Tm_Na_adjust / 1000.0) * math.log(Tm_Na_adjust / 1000.0, math.e)))
            d = 1.42 * pow(10, -5) * (
                    1.279 - 4.03 * pow(10, -3) * math.log(Tm_Na_adjust / 1000.0, math.e) - 8.03 * pow(10, -3) * pow(
                math.log(Tm_Na_adjust / 1000.0, math.e), 2))
            g = 8.31 * pow(10, -5) * (
                    0.486 - 0.258 * math.log(Tm_Na_adjust / 1000.0, math.e) + 5.25 * pow(10, -3) * pow(
                math.log(Tm_Na_adjust / 1000.0, math.e), 3))
        # Eq 16
        correction = a + (b * math.log(free_divalent, math.e))
        + GC_fraction(seq) * (c + (d * math.log(free_divalent, math.e)))
        + (1 / (2 * (len(seq) - 1))) * (e + (f * math.log(free_divalent, math.e))
                                        + g * (pow((math.log(free_divalent, math.e)), 2)))

    if symmetry(seq):
        # Equation A
        Tm = round(1 / ((1 / (delta_H / (delta_S + 1.9872 * math.log(primer_concentration / (1 * pow(10, 9)), math.e))))
                        + correction) - Kelvin, 2)
    else:
        # Equation B
        Tm = round(1 / ((1 / (delta_H / (delta_S + 1.9872 * math.log(primer_concentration / (4 * pow(10, 9)), math.e))))
                        + correction) - Kelvin, 2)
    return Tm


##############################################################################################


class NN_degenerate(object):
    def __init__(self, seq_file, primer_length=18, coverage=0.8, number_of_dege_bases=18, score_of_dege_bases=1000,
                 product_len=250, position="2,-1", variation=2, raw_entropy_threshold=3.6, distance=4, GC="0.4,0.6",
                 nproc=10, outfile=""):
        self.primer_length = primer_length  # primer length
        self.coverage = coverage  # min coverage
        self.number_of_dege_bases = number_of_dege_bases
        self.score_of_dege_bases = score_of_dege_bases
        self.product = product_len
        self.position = position  # gap position
        self.Y_strict, self.Y_strict_R = self.get_Y()
        self.variation = variation  # coverage of n-nt variation and max_gap_number
        self.distance = distance  # haripin
        self.GC = GC.split(",")
        self.nproc = nproc  # GC content
        self.seq_dict, self.total_sequence_number = self.parse_seq(seq_file)
        self.position_list = self.seq_attribute(self.seq_dict)
        self.start_position = self.position_list[0]
        self.stop_position = self.position_list[1]
        self.length = self.position_list[2]
        self.raw_entropy_threshold = raw_entropy_threshold
        self.entropy_threshold = self.entropy_threshold_adjust(self.length)
        self.outfile = outfile
        self.resQ = Manager().Queue()

    # expand degenerate primer into a list.
    @staticmethod
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
        return ["".join(i) for i in product(*seq)]

    ##################################################
    ################# pre_filter #####################
    ##################################################

    ################### hairpin ######################
    def hairpin_check(self, primer):
        n = 0
        distance = self.distance
        while n <= len(primer) - 5 - 5 - distance:
            kmer = self.degenerate_seq(primer[n:n + 5])
            left = self.degenerate_seq(primer[n + 5 + distance:])
            for k in kmer:
                for l in left:
                    if re.search(RC(k), l):
                        return True
            n += 1
        return False

    ################# GC content #####################
    def GC_fraction(self, sequence):
        sequence_expand = self.degenerate_seq(sequence)
        GC_list = []
        for seq in sequence_expand:
            GC_list.append(round((list(seq).count("G") + list(seq).count("C")) / len(list(seq)), 3))
        GC_average = round(mean(GC_list), 2)
        return GC_average

    ################# di_nucleotide #####################
    def di_nucleotide(self, primer):
        primers = self.degenerate_seq(primer)
        for m in primers:
            for n in di_nucleotides:
                if re.search(n, m):
                    return True
        return False

    ################## GC Clamp ######################
    def GC_clamp(self, primer, num=4, length=13):
        for i in range(num, (num + length)):
            s = primer[-i:]
            gc_fraction = self.GC_fraction(s)
            if gc_fraction > 0.6:
                return True
        return False

    ################# position of degenerate base #####################
    def dege_filter_in_term_N_bp(self, sequence):
        term = self.position
        if term == 0:
            term_base = ["A"]
        else:
            term_base = sequence[-term:]
        score = score_trans(term_base)
        if score > 1:
            return True
        else:
            return False

    # Import multi-alignment results and return a dict ==> {ID：sequence}
    def parse_seq(self, Input):
        seq_dict = defaultdict(str)
        with open(Input, "r") as f:
            for i in f:
                if i.startswith("#"):
                    pass
                else:
                    if i.startswith(">"):
                        i = i.strip().split(" ")
                        acc_id = i[0]
                    else:
                        # carefully !, make sure that Ns have been replaced!
                        sequence = re.sub("[^ACGTRYMKSWHBVD]", "-", i.strip().upper())
                        seq_dict[acc_id] += sequence
        return seq_dict, len(seq_dict)

    def current_end(self, primer, adaptor="", num=5, length=14):
        primer_extend = adaptor + primer
        end_seq = []
        for i in range(num, (num + length)):
            s = primer_extend[-i:]
            if s:
                end_seq.extend(self.degenerate_seq(s))
        return end_seq

    def deltaG(self, sequence):
        Delta_G_list = []
        Na = 50
        for seq in self.degenerate_seq(sequence):
            Delta_G = 0
            for n in range(len(seq) - 1):
                base_i, base_j = base2bit[seq[n + 1]], base2bit[seq[n]]
                Delta_G += freedom_of_H_37_table[base_i][base_j] * H_bonds_number[base_i][base_j] + \
                           penalty_of_H_37_table[base_i][base_j]
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

    def dimer_check(self, primer):
        current_end = self.current_end(primer)
        current_end_sort = sorted(current_end, key=lambda i: len(i), reverse=True)
        for end in current_end_sort:
            for p in self.degenerate_seq(primer):
                idx = p.find(RC(end))
                if idx >= 0:
                    end_length = len(end)
                    end_GC = end.count("G") + end.count("C")
                    end_d1 = 0
                    end_d2 = len(p) - len(end) - idx
                    Loss = Penalty_points(
                        end_length, end_GC, end_d1, end_d2)
                    delta_G = self.deltaG(end)
                    if Loss >= 3 or (delta_G < -5 and (end_d1 == end_d2)):
                        return True
        return False

    ####################################################################
    ##### pre-filter by GC content / di-nucleotide / hairpin ###########
    def primer_pre_filter(self, primer):
        information = []
        min_GC, max_GC = self.GC
        primer_GC_content = self.GC_fraction(primer)
        if not float(min_GC) <= primer_GC_content <= float(max_GC):
            information.append("GC_out_of_range (" + str(primer_GC_content) + ")")
        if self.di_nucleotide(primer):
            information.append("di_nucleotide")
        if self.hairpin_check(primer):
            information.append("hairpin")

        if len(information) == 0:
            return primer_GC_content
        else:
            return '|'.join(information)

    ####################################################################
    # if full degenerate primer is ok, we don't need to continue NN-array
    def pre_degenerate_primer_check(self, primer):
        primer_degeneracy = score_trans(primer)
        primer_dege_number = dege_number(primer)
        if primer_degeneracy < self.score_of_dege_bases and primer_dege_number < self.number_of_dege_bases:
            return True
        else:
            return False

    def full_degenerate_primer(self, freq_matrix):
        # degenerate transformation in each position
        max_dege_primers = ''
        for col in freq_matrix.columns.values:
            tmp = freq_matrix[freq_matrix[col] > 0].index.values.tolist()
            max_dege_primers += trans_score_table[round(sum([score_table[x] for x in tmp]), 2)]
        return max_dege_primers

    def state_matrix(self, primers_db):
        pieces = []
        for col in primers_db.columns.values:
            tmp_series = primers_db[col].value_counts()  # (normalize=True)
            tmp_series.name = col
            pieces.append(tmp_series)
        nodes = pd.concat(pieces, axis=1)
        nodes.fillna(0, inplace=True)
        nodes = nodes.sort_index(ascending=True)
        nodes = nodes.astype(int)
        row_names = nodes.index.values.tolist()
        if "-" in row_names:
            nodes.drop("-", inplace=True, axis=0)
        return nodes

    def di_matrix(self, primers):
        primers_trans = []
        for i in primers.keys():
            slice = []
            for j in range(len(i) - 1):
                slice.append(i[j:j + 2])
            primers_trans.extend(repeat(slice, primers[i]))
        return pd.DataFrame(primers_trans)

    def trans_matrix(self, primers):
        primers_di_db = self.di_matrix(primers)
        pieces = []
        for col in primers_di_db.columns.values:
            tmp_list = []
            for i in di_bases:
                # row: A, T, C ,G; column: A, T, C, G
                number = list(primers_di_db[col]).count(i)
                tmp_list.append(number)
            pieces.append(tmp_list)
        a, b = primers_di_db.shape
        trans = np.array(pieces).reshape(b, 4, 4)
        return trans

    def get_optimal_primer_by_viterbi(self, nodes, trans):
        nodes = np.array(nodes.T)
        seq_len, num_labels = len(nodes), len(trans[0])
        labels = np.arange(num_labels).reshape((1, -1))
        scores = nodes[0].reshape((-1, 1))
        primer_index = labels
        for t in range(1, seq_len):
            observe = nodes[t].reshape((1, -1))
            current_trans = trans[t - 1]
            M = scores + current_trans + observe
            scores = np.max(M, axis=0).reshape((-1, 1))
            idxs = np.argmax(M, axis=0)
            primer_index = np.concatenate([primer_index[:, idxs], labels], 0)
        best_primer_index = primer_index[:, scores.argmax()]
        return best_primer_index

    def get_optimal_primer_by_MM(self, cover_for_MM):
        sort_cover = sorted(cover_for_MM.items(), key=lambda x: x[1], reverse=True)
        L_seq = list(sort_cover[0][0])
        best_primer_index = [base2bit[x] for x in L_seq]
        # Return the maximum of an array or maximum along an axis. axis=0 代表列 , axis=1 代表行
        return best_primer_index

    def entropy(self, cover, cover_number, gap_sequence, gap_sequence_number):
        # cBit: entropy of cover sequences
        # tBit: entropy of total sequences
        cBit = 0
        tBit = 0
        for c in cover.keys():
            cBit += (cover[c] / cover_number) * math.log((cover[c] / cover_number), 2)
            tBit += (cover[c] / (cover_number + gap_sequence_number)) * \
                    math.log((cover[c] / (cover_number + gap_sequence_number)), 2)
        for t in gap_sequence.keys():
            tBit += (gap_sequence[t] / (cover_number + gap_sequence_number)) * \
                    math.log((gap_sequence[t] / (cover_number + gap_sequence_number)), 2)
        return round(-cBit, 2), round(-tBit, 2)

    # Sequence processing. Return a list contains sequence length, start and stop position of each sequence.
    def seq_attribute(self, Input_dict):
        start_dict = {}
        stop_dict = {}
        # pattern_start = re.compile('[A-Z]')
        # pattern_stop = re.compile("-*$")
        for acc_id in Input_dict.keys():
            # start_dict[acc_id] = pattern_start.search(Input_dict[acc_id]).span()[0]
            # stop_dict[acc_id] = pattern_stop.search(Input_dict[acc_id]).span()[0] - 1
            t_length = len(Input_dict[acc_id])
            start_dict[acc_id] = t_length - len(Input_dict[acc_id].lstrip("-"))
            stop_dict[acc_id] = len(Input_dict[acc_id].rstrip("-"))
            # start position should contain [coverage] sequences at least.
        start = np.quantile(np.array(list(start_dict.values())).reshape(1, -1), self.coverage, interpolation="higher")
        # for python 3.9.9
        # start = np.quantile(np.array(list(start_dict.values())).reshape(1, -1), self.coverage, method="higher")
        # stop position should contain [coverage] sequences at least.
        stop = np.quantile(np.array(list(stop_dict.values())).reshape(1, -1), self.coverage, interpolation="lower")
        # stop = np.quantile(np.array(list(stop_dict.values())).reshape(1, -1), self.coverage, method="lower")
        if stop - start < int(self.product):
            print("Error: max length of PCR product is shorter than the default min Product length with {} "
                  "coverage! Non candidate primers !!!".format(self.coverage))
            sys.exit(1)
        else:
            return [start, stop, stop - start]

    def entropy_threshold_adjust(self, length):
        if length < 5000:
            return self.raw_entropy_threshold
        else:
            if length < 10000:
                return self.raw_entropy_threshold * 0.95
            else:
                return self.raw_entropy_threshold * 0.9

    def get_primers(self, sequence_dict, primer_start):  # , primer_info, non_cov_primer_out
        # record sequence and acc id
        non_gap_seq_id = defaultdict(list)
        # record sequence (no gap) and number
        cover = defaultdict(int)
        cover_for_MM = defaultdict(int)
        # record total coverage sequence number
        cover_number = 0
        # record sequence (> variation gap) and number
        gap_sequence = defaultdict(int)
        gap_seq_id = defaultdict(list)
        # record total sequence (> variation gap) number
        gap_sequence_number = 0
        primers_db = []
        for seq_id in sequence_dict.keys():
            sequence = sequence_dict[seq_id][primer_start:primer_start + self.primer_length].upper()
            # replace "-" which in start or stop position with nucleotides
            if sequence == "-" * self.primer_length:
                pass
            else:
                if sequence.startswith("-"):
                    sequence_narrow = sequence.lstrip("-")
                    append_base_length = len(sequence) - len(sequence_narrow)
                    left_seq = sequence_dict[seq_id][0:primer_start].replace("-", "")
                    if len(left_seq) >= append_base_length:
                        sequence = left_seq[len(left_seq) - append_base_length:] + sequence_narrow
                if sequence.endswith("-"):
                    sequence_narrow = sequence.rstrip("-")
                    append_base_length = len(sequence) - len(sequence_narrow)
                    right_seq = sequence_dict[seq_id][primer_start + self.primer_length:].replace("-", "")
                    if len(right_seq) >= append_base_length:
                        sequence = sequence_narrow + right_seq[0:append_base_length]
            if len(sequence) < self.primer_length:
                append_base_length = self.primer_length - len(sequence)
                left_seq = sequence_dict[seq_id][0:primer_start].replace("-", "")
                if len(left_seq) >= append_base_length:
                    sequence = left_seq[len(left_seq) - append_base_length:] + sequence
            # gap number. number of gap > 2
            if list(sequence).count("-") > self.variation:
                gap_sequence[sequence] += 1
                gap_sequence_number += 1
                if round(gap_sequence_number / self.total_sequence_number, 2) >= (1 - self.coverage):
                    break
                else:
                    # record acc ID of gap sequences
                    expand_sequence = self.degenerate_seq(sequence)
                    for i in expand_sequence:
                        gap_seq_id[i].append(seq_id)
            # # accepted gap, number of gap <= variation
            else:
                expand_sequence = self.degenerate_seq(sequence)
                cover_number += 1
                for i in expand_sequence:
                    cover[i] += 1
                    primers_db.append(list(i))
                    # record acc ID of non gap sequences, which is potential mis-coverage
                    non_gap_seq_id[i].append(seq_id)
                    if re.search("-", i):
                        pass
                    else:
                        cover_for_MM[i] += 1
        # number of sequences with too many gaps greater than (1 - self.coverage)
        if round(gap_sequence_number / self.total_sequence_number, 2) >= (1 - self.coverage):
            # print("Gap fail")
            self.resQ.put(None)
        elif len(cover) < 1:
            self.resQ.put(None)
            # print("Cover fail")
        else:
            # cBit: entropy of cover sequences
            # tBit: entropy of total sequences
            cBit, tBit = self.entropy(cover, cover_number, gap_sequence, gap_sequence_number)
            if tBit > self.entropy_threshold:
                # print("Entropy fail")
                # This window is not a conserved region, and not proper to design primers
                self.resQ.put(None)
            else:
                primers_db = pd.DataFrame(primers_db)
                # frequency matrix
                freq_matrix = self.state_matrix(primers_db)
                # print(freq_matrix)
                colSum = np.sum(freq_matrix, axis=0)
                a, b = freq_matrix.shape
                # a < 4 means base composition of this region is less than 4 (GC bias).
                # It's not a proper region for primer design.
                if a < 4:
                    self.resQ.put(None)
                elif (colSum == 0).any():
                    # print(colSum)  # if 0 in array; pass
                    self.resQ.put(None)
                else:
                    gap_seq_id_info = [primer_start, gap_seq_id]
                    mismatch_coverage, non_cov_primer_info = \
                        self.degenerate_by_NN_algorithm(primer_start, freq_matrix, cover, non_gap_seq_id,
                                                        cover_for_MM, cover_number, cBit, tBit)
                    # self.resQ.put([mismatch_coverage, non_cov_primer_info, gap_seq_id_info])
                    # F, R = mismatch_coverage[1][6], mismatch_coverage[1][7]
                    sequence = mismatch_coverage[1][2]
                    if self.dimer_check(sequence):
                        # print("Dimer fail")
                        self.resQ.put(None)
                    else:
                        self.resQ.put([mismatch_coverage, non_cov_primer_info, gap_seq_id_info])
                    # if F < cover_number * 0.5 or R < cover_number * 0.5:
                    #     self.resQ.put(None)
                    # else:
                    #     self.resQ.put([mismatch_coverage, non_cov_primer_info, gap_seq_id_info])

    def degenerate_by_NN_algorithm(self, primer_start, freq_matrix, cover, non_gap_seq_id, cover_for_MM,
                                   cover_number, cBit, tBit):
        # full degenerate primer
        full_degenerate_primer = self.full_degenerate_primer(freq_matrix)
        # unique covered primers, which is used to calculate coverage and
        # mis-coverage in the following step.
        cover_primer_set = set(cover.keys())
        # if full_degenerate_primer is ok, then return full_degenerate_primer
        # mismatch_coverage, non_cov_primer_info = {}, {}
        F_non_cover, R_non_cover = {}, {}
        ######################################################################################################
        ############ need prone. not all primers is proper for primer-F or primer-R ##########################
        # If a primer is located in the start region, there is no need to calculate its coverage for primer-R#
        ## here is a suggestion. we can assert candidate primer as primer-F or primer-R by primer attribute ##
        ######################################################################################################
        NN_matrix = self.trans_matrix(cover)
        if len(cover_for_MM) != 0:
            optimal_primer_index_NM = self.get_optimal_primer_by_viterbi(freq_matrix, NN_matrix)
            optimal_primer_index_MM = self.get_optimal_primer_by_MM(cover_for_MM)
            # print(optimal_primer_index_NM.tolist()) # array
            # print(optimal_primer_index_MM) # list
            #  if (optimal_primer_index_NM == optimal_primer_index_MM).all():
            if optimal_primer_index_NM.tolist() == optimal_primer_index_MM:
                optimal_primer_index = optimal_primer_index_NM
                row_names = np.array(freq_matrix.index.values).reshape(1, -1)
                # build a list to store init base information in each position.
                optimal_primer_list = row_names[:, optimal_primer_index][0].tolist()
                # initiation coverage (optimal primer, used as base coverage)
                optimal_coverage_init = cover["".join(optimal_primer_list)]
                optimal_primer_current, F_mis_cover, R_mis_cover, information, F_non_cover, R_non_cover = \
                    self.coverage_stast(cover, optimal_primer_index, NN_matrix, optimal_coverage_init, cover_number,
                                        optimal_primer_list, cover_primer_set, non_gap_seq_id, F_non_cover,
                                        R_non_cover)
                # print(F_mis_cover)
                # print(R_mis_cover)
            else:
                F_non_cover_NM, R_non_cover_NM, F_non_cover_MM, R_non_cover_MM = {}, {}, {}, {}
                row_names = np.array(freq_matrix.index.values).reshape(1, -1)
                # build a list to store init base information in each position.
                optimal_primer_list_NM = row_names[:, optimal_primer_index_NM][0].tolist()
                # initiation coverage (optimal primer, used as base coverage)
                optimal_coverage_init_NM = cover["".join(optimal_primer_list_NM)]
                NN_matrix_NM = NN_matrix.copy()
                optimal_primer_current_NM, F_mis_cover_NM, R_mis_cover_NM, information_NM, F_non_cover_NM, \
                R_non_cover_NM = self.coverage_stast(cover, optimal_primer_index_NM, NN_matrix_NM,
                                                     optimal_coverage_init_NM, cover_number, optimal_primer_list_NM,
                                                     cover_primer_set, non_gap_seq_id, F_non_cover_NM,
                                                     R_non_cover_NM)
                optimal_primer_list_MM = row_names[:, optimal_primer_index_MM][0].tolist()
                # initiation coverage (optimal primer, used as base coverage)
                optimal_coverage_init_MM = cover["".join(optimal_primer_list_MM)]
                NN_matrix_MM = NN_matrix.copy()
                optimal_primer_current_MM, F_mis_cover_MM, R_mis_cover_MM, information_MM, F_non_cover_MM, \
                R_non_cover_MM = self.coverage_stast(cover, optimal_primer_index_MM, NN_matrix_MM,
                                                     optimal_coverage_init_MM, cover_number,
                                                     optimal_primer_list_MM, cover_primer_set, non_gap_seq_id,
                                                     F_non_cover_MM, R_non_cover_MM)
                if (F_mis_cover_NM + R_mis_cover_NM) > (F_mis_cover_MM + R_mis_cover_MM):
                    optimal_primer_current, F_mis_cover, R_mis_cover, information, optimal_coverage_init, \
                    F_non_cover, R_non_cover, NN_matrix = optimal_primer_current_NM, F_mis_cover_NM, \
                                                          R_mis_cover_NM, information_NM, optimal_coverage_init_NM, \
                                                          F_non_cover_NM, R_non_cover_NM, NN_matrix_NM
                else:
                    optimal_primer_current, F_mis_cover, R_mis_cover, information, optimal_coverage_init, \
                    F_non_cover, R_non_cover, NN_matrix = optimal_primer_current_MM, F_mis_cover_MM, \
                                                          R_mis_cover_MM, information_MM, optimal_coverage_init_MM, \
                                                          F_non_cover_MM, R_non_cover_MM, NN_matrix_MM
                # print(F_mis_cover)
                # print(R_mis_cover)
        else:
            optimal_primer_index_NM = self.get_optimal_primer_by_viterbi(freq_matrix, NN_matrix)
            F_non_cover_NM, R_non_cover_NM, F_non_cover_MM, R_non_cover_MM = {}, {}, {}, {}
            row_names = np.array(freq_matrix.index.values).reshape(1, -1)
            # build a list to store init base information in each position.
            optimal_primer_list_NM = row_names[:, optimal_primer_index_NM][0].tolist()
            # initiation coverage (optimal primer, used as base coverage)
            optimal_coverage_init_NM = cover["".join(optimal_primer_list_NM)]
            NN_matrix_NM = NN_matrix.copy()
            optimal_primer_current_NM, F_mis_cover_NM, R_mis_cover_NM, information_NM, F_non_cover_NM, \
            R_non_cover_NM = self.coverage_stast(cover, optimal_primer_index_NM, NN_matrix_NM,
                                                 optimal_coverage_init_NM, cover_number, optimal_primer_list_NM,
                                                 cover_primer_set, non_gap_seq_id, F_non_cover_NM, R_non_cover_NM)
            optimal_primer_current, F_mis_cover, R_mis_cover, information, optimal_coverage_init, F_non_cover, \
            R_non_cover, NN_matrix = optimal_primer_current_NM, F_mis_cover_NM, R_mis_cover_NM, information_NM, \
                                     optimal_coverage_init_NM, F_non_cover_NM, R_non_cover_NM, NN_matrix_NM
            # print(F_mis_cover)
            # print(R_mis_cover)
        nonsense_primer_number = len(set(self.degenerate_seq(optimal_primer_current)) - set(cover.keys()))
        primer_degenerate_number = dege_number(optimal_primer_current)
        Tm, coverage = [], []
        for seq in self.degenerate_seq(optimal_primer_current):
            Tm.append(Calc_Tm_v2(seq))
            coverage.append(cover[seq])
        Tm_average = round(mean(Tm), 2)
        perfect_coverage = sum(coverage)
        out_mismatch_coverage = [primer_start, [cBit, tBit, optimal_primer_current, primer_degenerate_number,
                                                nonsense_primer_number, perfect_coverage, F_mis_cover,
                                                R_mis_cover, Tm_average, information]]
        non_cov_primer_info = [primer_start, [F_non_cover, R_non_cover]]
        return out_mismatch_coverage, non_cov_primer_info

    def coverage_stast(self, cover, optimal_primer_index, NN_matrix, optimal_coverage_init, cover_number,
                       optimal_primer_list, cover_primer_set, non_gap_seq_id, F_non_cover, R_non_cover):
        # if the coverage is too low, is it necessary to refine?
        # mis-coverage as threshold? if mis-coverage reached to 100% but degeneracy is still very low,
        optimal_NN_index = []
        optimal_NN_coverage = []
        for idx in range(len(optimal_primer_index) - 1):
            # NN index
            optimal_NN_index.append([optimal_primer_index[idx], optimal_primer_index[idx + 1]])
            # NN coverage
            # Is the minimum number in NN coverage = optimal_primer_coverage ? No!
            optimal_NN_coverage.append(
                NN_matrix[idx, optimal_primer_index[idx], optimal_primer_index[idx + 1]])
        # mis-coverage initialization
        F_mis_cover_cover, F_non_cover_in_cover, R_mis_cover_cover, R_non_cover_in_cover = \
            self.mis_primer_check(cover_primer_set, ''.join(optimal_primer_list), cover,
                                  non_gap_seq_id)
        # print(optimal_coverage_init + F_mis_cover_cover)
        # print(optimal_coverage_init + R_mis_cover_cover)
        # print(cover_number)
        # print(optimal_primer_list)
        if optimal_coverage_init + F_mis_cover_cover < cover_number or \
                optimal_coverage_init + R_mis_cover_cover < cover_number:
            while optimal_coverage_init + F_mis_cover_cover < cover_number or \
                    optimal_coverage_init + R_mis_cover_cover < cover_number:
                # optimal_primer_update, coverage_update, NN_coverage_update,
                # NN array_update, degeneracy_update, degenerate_update
                optimal_primer_list, optimal_coverage_init, optimal_NN_coverage_update, \
                NN_matrix, degeneracy, number_of_degenerate = \
                    self.refine_by_NN_array(optimal_primer_list, optimal_coverage_init, cover, optimal_NN_index,
                                            optimal_NN_coverage, NN_matrix)
                F_mis_cover_cover, F_non_cover_in_cover, R_mis_cover_cover, R_non_cover_in_cover = \
                    self.mis_primer_check(cover_primer_set, ''.join(optimal_primer_list), cover,
                                          non_gap_seq_id)
                # If there is no increase in NN_coverage,
                # it suggests the presence of bugs or a mismatch in continuous positions.
                # Is this step necessary? or shall we use DegePrime method? or shall we use machine learning?
                if max(F_mis_cover_cover, R_mis_cover_cover) == cover_number:
                    break
                elif optimal_NN_coverage_update == optimal_NN_coverage:
                    break
                # If the degeneracy exceeds the threshold, the loop will break.
                elif 2 * degeneracy > self.score_of_dege_bases or 3 * degeneracy / 2 > self.score_of_dege_bases \
                        or number_of_degenerate == self.number_of_dege_bases:
                    break
                else:
                    optimal_NN_coverage = optimal_NN_coverage_update
        # If the primer coverage does not increase after degeneration,
        # the process will backtrack and assess the original optimal primer.
        # print(optimal_primer_list)
        optimal_primer_current = ''.join(optimal_primer_list)
        information = self.primer_pre_filter(optimal_primer_current)
        # F_mis_cover_cover, F_non_cover_in_cover, R_mis_cover_cover, R_non_cover_in_cover = \
        #     self.mis_primer_check(cover_primer_set, optimal_primer_current, cover,
        #                           non_gap_seq_id)
        F_non_cover.update(F_non_cover_in_cover)
        R_non_cover.update(R_non_cover_in_cover)
        F_mis_cover = optimal_coverage_init + F_mis_cover_cover
        R_mis_cover = optimal_coverage_init + R_mis_cover_cover
        # print(F_mis_cover)
        return optimal_primer_current, F_mis_cover, R_mis_cover, information, F_non_cover, R_non_cover

    def refine_by_NN_array(self, optimal_primer_list, optimal_coverage_init, cover,
                           optimal_NN_index, optimal_NN_coverage, NN_array):
        # use minimum index of optimal_NN_coverage as the position to refine
        refine_index = np.where(optimal_NN_coverage == np.min(optimal_NN_coverage))[0]  # np.where[0] is a list
        # build dict to record coverage and NN array
        primer_update_list, coverage_update_list, NN_array_update_list, NN_coverage_update = [], [], [], []
        for i in refine_index:
            optimal_NN_coverage_tmp = optimal_NN_coverage.copy()
            NN_array_tmp = NN_array.copy()
            optimal_list = optimal_primer_list.copy()
            # initiation score
            # initiation coverage
            coverage_renew = optimal_coverage_init
            if i == 0:
                # two position need refine
                # position 0 and 1
                # decide which position to choose
                row = optimal_NN_index[i][0]
                column = optimal_NN_index[i][1]
                if len(np.where(NN_array_tmp[0, :, column] > 0)[0]) > 1:
                    init_score = score_table[optimal_list[i]]
                    refine_column = NN_array_tmp[i, :, column]
                    refine_row_arg_sort = np.argsort(refine_column, axis=0)[::-1]
                    new_primer = optimal_list
                    # print(row, refine_column)
                    for idx in refine_row_arg_sort:
                        # init refine,  We must ensure that there are no double counting.
                        # position 0.
                        if idx != row:
                            init_score += score_table[bases[idx]]
                            new_primer[i] = bases[idx]
                            # Calculate coverage after refine
                            for new_primer_update in self.degenerate_seq("".join(new_primer)):
                                if new_primer_update in cover.keys():
                                    coverage_renew += cover["".join(new_primer_update)]
                            new_primer[i] = trans_score_table[round(init_score, 2)]
                            # reset NN_array. row names will update after reset.
                            NN_array_tmp[i, row, :] += NN_array_tmp[i, idx, :]
                            NN_array_tmp[i, idx, :] -= NN_array_tmp[i, idx, :]
                            optimal_NN_coverage_tmp[i] = NN_array_tmp[i, row, column]
                            break
                        # primer update
                    optimal_list_update = optimal_list
                    optimal_list_update[i] = trans_score_table[round(init_score, 2)]
                # position 1
                elif len(np.where(NN_array_tmp[0, row, :] > 0)[0]) > 1:
                    init_score = score_table[optimal_list[i + 1]]
                    next_row = optimal_NN_index[i + 1][0]
                    next_column = optimal_NN_index[i + 1][1]
                    # concat row of layer i and column of layer i+1
                    refine_row = NN_array_tmp[i, row, :].reshape(1, -1)
                    refine_column = NN_array_tmp[i + 1, :, next_column].reshape(1, -1)
                    refine = np.concatenate([refine_row, refine_column], 0)
                    refine_min = np.min(refine, axis=0)
                    refine_row_arg_sort = np.argsort(refine_min, axis=0)[::-1]
                    # Return the minimum of an array or maximum along an axis. axis=0: column , axis=1: row
                    new_primer = optimal_list
                    if len(np.where(refine_min > 0)[0]) > 1:
                        for idx in refine_row_arg_sort:
                            # We must ensure that there are no double counting.
                            # position 1.
                            if idx != column:
                                init_score += score_table[bases[idx]]
                                # Calculate coverage after refine
                                new_primer[i + 1] = bases[idx]
                                for new_primer_update in self.degenerate_seq("".join(new_primer)):
                                    if new_primer_update in cover.keys():
                                        coverage_renew += cover["".join(new_primer_update)]
                                new_primer[i + 1] = trans_score_table[round(init_score, 2)]
                                # reset NN_array. column + (column idx) of layer i and row + (row idx) of layer i+1.
                                NN_array_tmp[i, :, column] += NN_array_tmp[i, :, idx]
                                NN_array_tmp[i, :, idx] -= NN_array_tmp[i, :, idx]
                                NN_array_tmp[i + 1, next_row, :] += NN_array_tmp[i + 1, idx, :]
                                NN_array_tmp[i + 1, idx, :] -= NN_array_tmp[i + 1, idx, :]
                                optimal_NN_coverage_tmp[i] = NN_array_tmp[i, row, column]
                                optimal_NN_coverage_tmp[i + 1] = NN_array_tmp[i + 1, next_row, next_column]
                                break
                    # primer update
                    optimal_list_update = optimal_list
                    optimal_list_update[i + 1] = trans_score_table[round(init_score, 2)]
                else:
                    optimal_list_update = optimal_list
            elif i == len(optimal_NN_index) - 1:
                init_score = score_table[optimal_list[i + 1]]
                row = optimal_NN_index[i][0]
                column = optimal_NN_index[i][1]
                refine_row = NN_array_tmp[i, row, :]
                refine_row_arg_sort = np.argsort(refine_row, axis=0)[::-1]
                # If number of refine_row > 1, then the current position need to refine.
                if len(np.where(refine_row > 0)[0]) > 1:
                    new_primer = optimal_list
                    for idx in refine_row_arg_sort:
                        # We must ensure that there are no double counting.
                        # position -1.
                        if idx != column:
                            init_score += score_table[bases[idx]]
                            # Calculate coverage after refine
                            new_primer[i + 1] = bases[idx]
                            for new_primer_update in self.degenerate_seq("".join(new_primer)):
                                if new_primer_update in cover.keys():
                                    coverage_renew += cover["".join(new_primer_update)]
                            new_primer[i + 1] = trans_score_table[round(init_score, 2)]
                            # reset NN_array. column names will update after reset.
                            NN_array_tmp[i, :, column] += NN_array_tmp[i, :, idx]
                            NN_array_tmp[i, :, idx] -= NN_array_tmp[i, :, idx]
                            optimal_NN_coverage_tmp[i] = NN_array_tmp[i, row, column]
                            break
                # primer update
                optimal_list_update = optimal_list
                optimal_list_update[i + 1] = trans_score_table[round(init_score, 2)]
            else:
                init_score = score_table[optimal_list[i + 1]]
                row = optimal_NN_index[i][0]
                column = optimal_NN_index[i][1]
                next_row = optimal_NN_index[i + 1][0]
                next_column = optimal_NN_index[i + 1][1]
                # concat row of layer i and column of layer i+1
                refine_row = NN_array_tmp[i, row, :].reshape(1, -1)
                refine_column = NN_array_tmp[i + 1, :, next_column].reshape(1, -1)
                refine = np.concatenate([refine_row, refine_column], 0)
                refine_min = np.min(refine, axis=0)
                # Return the minimum of an array or maximum along an axis. axis=0: column , axis=1: row
                refine_min_arg_sort = np.argsort(refine_min, axis=0)[::-1]
                if len(np.where(refine_min > 0)[0]) > 1:
                    new_primer = optimal_list
                    # for idx in np.where(refine_min_sort > 0)[0]:
                    for idx in refine_min_arg_sort:
                        # We must ensure that there are no double counting.
                        # position i+1.
                        if idx != column:
                            # or if idx != next_row
                            # init trans score update
                            init_score += score_table[bases[idx]]
                            # Calculate coverage after refine
                            new_primer[i + 1] = bases[idx]
                            for new_primer_update in self.degenerate_seq("".join(new_primer)):
                                if new_primer_update in cover.keys():
                                    coverage_renew += cover["".join(new_primer_update)]
                            new_primer[i + 1] = trans_score_table[round(init_score, 2)]
                            # reset NN_array. column + (column idx) of layer i and row + (row idx) of layer i+1.
                            NN_array_tmp[i, :, column] += NN_array_tmp[i, :, idx]
                            NN_array_tmp[i, :, idx] -= NN_array_tmp[i, :, idx]
                            NN_array_tmp[i + 1, next_row, :] += NN_array_tmp[i + 1, idx, :]
                            NN_array_tmp[i + 1, idx, :] -= NN_array_tmp[i + 1, idx, :]
                            optimal_NN_coverage_tmp[i] = NN_array_tmp[i, row, column]
                            optimal_NN_coverage_tmp[i + 1] = NN_array_tmp[i + 1, next_row, next_column]

                            break
                # primer update
                optimal_list_update = optimal_list
                optimal_list_update[i + 1] = trans_score_table[round(init_score, 2)]
            # primer_update = "".join(primer_list_update)
            primer_update_list.append(optimal_list_update)
            NN_coverage_update.append(optimal_NN_coverage_tmp)
            # current_primers_set = set(self.degenerate_seq(primer_update))
            # coverage of update primers
            coverage_update_list.append(coverage_renew)
            # new NN_array
            NN_array_update_list.append(NN_array_tmp)
        optimal_idx = coverage_update_list.index(max(coverage_update_list))
        degeneracy_update = score_trans(primer_update_list[optimal_idx])
        degenerate_number_update = sum([math.floor(score_table[x]) > 1 for x in primer_update_list[optimal_idx]])
        # optimal_primer_update, coverage_update,
        # NN_coverage_update, NN array_update,
        # degeneracy_update, degenerate_update
        return primer_update_list[optimal_idx], coverage_update_list[optimal_idx], \
               NN_coverage_update[optimal_idx], NN_array_update_list[optimal_idx], \
               degeneracy_update, degenerate_number_update

    def get_Y(self):
        Y_strict, Y_strict_R = [], []
        for y in self.position.split(","):
            y_index = int(y.strip())
            if y_index > 0:
                Y_strict.append(y_index)
                Y_strict_R.append(self.primer_length - y_index)
            else:
                Y_strict.append(self.primer_length + y_index + 1)
                Y_strict_R.append(-y_index + 1)
        return set(Y_strict), set(Y_strict_R)

    def mis_primer_check(self, all_primers, optimal_primer, cover, non_gap_seq_id):
        # uncoverage sequence in cover dict
        optimal_primer_set = set(self.degenerate_seq(optimal_primer))
        uncover_primer_set = all_primers - optimal_primer_set
        F_non_cover, R_non_cover = {}, {}
        F_mis_cover, R_mis_cover = 0, 0
        for uncover_primer in uncover_primer_set:
            Y_dist = Y_distance(optimal_primer, uncover_primer)
            # print(uncover_primer)
            # print(Y_dist)
            # print(set(Y_dist))
            if len(Y_dist) > self.variation:
                # record sequence and acc_ID which will never mis-coverage. too many mismatch!
                F_non_cover[uncover_primer] = non_gap_seq_id[uncover_primer]
                R_non_cover[uncover_primer] = non_gap_seq_id[uncover_primer]
            # if len(Y_dist) <= self.variation:
            else:
                if len(set(Y_dist).intersection(self.Y_strict)) > 0:
                    F_non_cover[uncover_primer] = non_gap_seq_id[uncover_primer]
                else:
                    F_mis_cover += cover[uncover_primer]
                if len(set(Y_dist).intersection(self.Y_strict_R)) > 0:
                    R_non_cover[uncover_primer] = non_gap_seq_id[uncover_primer]
                else:
                    R_mis_cover += cover[uncover_primer]
        # print(optimal_primer)
        # print(F_mis_cover)
        return F_mis_cover, F_non_cover, R_mis_cover, R_non_cover

    ################# get_primers #####################
    def run(self):
        p = ProcessPoolExecutor(self.nproc)  #
        sequence_dict = self.seq_dict
        start_primer = self.start_position
        stop_primer = self.stop_position
        # primer_info = Manager().list()
        # non_cov_primer_out = Manager().list()
        # for position in range(1245,  stop_primer - self.primer_length):
        for position in range(start_primer, stop_primer - self.primer_length):
            # print(position)
            p.submit(self.get_primers(sequence_dict, position))  # , primer_info, non_cov_primer_out
            # This will submit all tasks to one place without blocking, and then each
            # thread in the thread pool will fetch tasks.
        n = 0
        candidate_list, non_cov_primer_out, gap_seq_id_out = [], [], []
        with open(self.outfile, "w") as fo:
            headers = ["Position", "Entropy of cover (bit)", "Entropy of total (bit)", "Optimal_primer",
                       "primer_degenerate_number",
                       "nonsense_primer_number", "Optimal_coverage", "Mis-F-coverage", "Mis-R-coverage", "Tm",
                       "Information"]
            fo.write("\t".join(map(str, headers)) + "\n")
            while n < stop_primer - start_primer - self.primer_length:
                res = self.resQ.get()
                # The get method can read and delete an element from the queue. Similarly, the get method has two
                # optional parameters: blocked and timeout. If blocked is true (the default value) and timeout is
                # positive, no element is retrieved during the waiting time, and a Queue is thrown Empty exception.
                # If blocked is false, there are two cases. If a value of Queue is available, return the value
                # immediately. Otherwise, if the queue is empty, throw a Queue.Empty exception immediately.
                if res is None:
                    n += 1
                    continue
                candidate_list.append(res[0])
                non_cov_primer_out.append(res[1])
                gap_seq_id_out.append(res[2])
                n += 1
            sorted_candidate_dict = dict(sorted(dict(candidate_list).items(), key=lambda x: x[0], reverse=False))
            for position in sorted_candidate_dict.keys():
                fo.write(str(position) + "\t" + "\t".join(map(str, sorted_candidate_dict[position])) + "\n")
            fo.close()
            with open(self.outfile + '.non_coverage_seq_id_json', "w") as fj:
                json.dump(dict(non_cov_primer_out), fj, indent=4)
            fj.close()
            with open(self.outfile + '.gap_seq_id_json', "w") as fg:
                json.dump(dict(gap_seq_id_out), fg, indent=4)
            fg.close()
            # get results before shutdown. Synchronous call mode: call, wait for the return value, decouple,
            # but slow.
        p.shutdown()


def main():
    args = parseArg()
    NN_APP = NN_degenerate(seq_file=args.input, primer_length=args.plen, coverage=args.fraction,
                           number_of_dege_bases=args.dnum, score_of_dege_bases=args.degeneracy,
                           raw_entropy_threshold=args.entropy, product_len=args.size, position=args.coordinate,
                           variation=args.variation, distance=args.away, GC=args.gc,
                           nproc=args.proc, outfile=args.out)
    NN_APP.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
