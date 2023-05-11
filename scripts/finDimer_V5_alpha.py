#!/bin/python

__date__ = "2022-10-8"
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
import os
import argparse
import time
from math import log10
from itertools import product
# from multiprocessing import Manager
from collections import defaultdict
# from concurrent.futures import ProcessPoolExecutor
import asyncio
import uvloop

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


def reversecomplement(seq):
    return seq.translate(TRANS)[::-1]


def Penalty_points(length, GC, d1, d2):
    # return log10((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)))
    return log10((2 ** length * 2 ** GC) / ((2 ** d1 - 0.9) * (2 ** d2 - 0.9)))

def parseArg():
    parser = argparse.ArgumentParser(
        description="For primer dimer check")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="input fasta primer file", metavar="<file>")
    parser.add_argument("-n", "--num", type=int, default=5,
                        help='number of async process, 5 by default', metavar="<int>")
    parser.add_argument("-t", "--threshold", type=float, default=3.96,
                        help='threshold of loss function. Default: 3.96', metavar="<int>")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help='output file', metavar="<file>")
    return parser.parse_args()


TRANS_c = str.maketrans("ATCG", "TAGC")


def complement(seq):
    return seq.translate(TRANS_c)[::-1]


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


class Dimer(object):

    def __init__(self, primer_file="", threshold=3.96, output_file='', nproc=10):
        self.nproc = nproc
        self.output_file = os.path.abspath(output_file)
        self.primers_file = primer_file
        self.threshold = threshold
        self.primers = self.parse_primers()
        self.sema = asyncio.Semaphore(self.nproc)
        self.primers_list = list(self.primers.keys())

    def parse_primers(self):
        primer_dict = defaultdict(str)
        with open(self.primers_file, "r") as f:
            for i in f:
                if i.startswith(">"):
                    name = i.strip()
                else:
                    primer_dict[i.strip()] = name
        return primer_dict

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

    def dimer_check(self, position):
        current_end_set = self.current_end(self.primers_list[position])
        current_end_list = sorted(list(current_end_set), key=lambda i: len(i), reverse=True)
        dimer_result = []
        dimer = False
        for ps in self.primers_list[position:]:
            for end in current_end_list:
                for p in self.degenerate_seq(ps):
                    idx = p.find(reversecomplement(end))
                    if idx >= 0:
                        end_length = len(end)
                        end_GC = end.count("G") + end.count("C")
                        end_d1 = 0
                        end_d2 = len(p) - len(end) - idx
                        Loss = Penalty_points(
                            end_length, end_GC, end_d1, end_d2)
                        delta_G = self.deltaG(end)
                        if Loss >= self.threshold or (delta_G < -5 and (end_d1 == end_d2)):
                            line = (self.primers[self.primers_list[position]], self.primers_list[position],
                                    end, delta_G, end_length, end_d1,
                                    end_GC, self.primers[ps],
                                    ps, end_d2, Loss
                                    )
                            dimer_result.append(line)
                            dimer = True
                            break
                if dimer:
                    dimer = False
                    break
        if dimer_result:
            # print(dimer_result)
            return dimer_result

    async def run(self, position, sema):
        async with sema:
            # print("This is {} subprocess".format(position))
            return self.dimer_check(position)

    def find_dimers(self):
        tasks = [self.run(position, self.sema) for position in range(len(self.primers_list))]
        asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())
        loop = asyncio.get_event_loop()
        results_tmp = loop.run_until_complete(asyncio.gather(*tasks))
        filtered_results = list(filter(lambda x: x is not None, results_tmp))
        results_item = [item for sublist in filtered_results for item in sublist]
        return results_item

    def write_results(self, results, output_file):
        primer_id_sum = defaultdict(int)
        dimer_primer_id_sum = defaultdict(int)
        with open(output_file, "w") as fo:
            headers = ["Primer_ID", "Primer seq", "Primer end", "Delta G", "Primer end length", "End (distance 1)",
                       "End (GC)", "Dimer-primer_ID", "Dimer-primer seq", "End (distance 2)", "Loss"]
            fo.write("\t".join(headers) + "\n")
            for result in results:
                primer_id_sum[result[0]] += 1
                dimer_primer_id_sum[result[7]] += 1
                fo.write("\t".join(map(str, result)) + "\n")

        with open(output_file + ".dimer_num", "w") as fo:
            fo.write("SeqName\tPrimer_ID\tDimer-primer_ID\tRowSum\n")
            for k in primer_id_sum.keys():
                p_id = primer_id_sum[k]
                d_id = dimer_primer_id_sum[k]
                RowSum = p_id + d_id
                fo.write("\t".join(map(str, [k, p_id, d_id, RowSum])) + "\n")

    def main(self):
        results = self.find_dimers()
        self.write_results(results, self.output_file)

def main():
    args = parseArg()
    dimer_app = Dimer(primer_file=args.input, threshold=args.threshold, nproc=args.num, output_file=args.output)
    dimer_app.main()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
