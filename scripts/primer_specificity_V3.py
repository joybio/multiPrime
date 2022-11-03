#!/bin/python
# coding:utf-8
"""
Output off-target PCR results. 
"""

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

import re
import sys
from collections import defaultdict
import os
import multiprocessing
from itertools import product
from multiprocessing import Process
import time
import numpy as np
import pandas as pd
from optparse import OptionParser
from pathlib import Path
import math
from math import log10
from functools import reduce
from operator import mul  #
from bisect import bisect_left
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor


# Path(path).parent, Path(path).name, Path(path).suffix, Path(path).stem, Path(path).iterdir(), Path(path).joinpath()
# Path(path).is_absolute(), Path(path).is_dir(), Path(path).is_file(), Path(path).exists(), Path(path).with_suffix


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -r [bowtie index] -l [150,2000] -p [10]-o [output]')

    parser.add_option('-i', '--input',
                      dest='input_file',
                      help='input file: primer.fa.')

    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='reference file: bowtie index.')

    parser.add_option('-l', '--len',
                      dest='len',
                      default=9,
                      type="int",
                      help='Length of primer term, which is used for blast. Default:9')

    parser.add_option('-s', '--s',
                      dest='size',
                      default="150,2000",
                      help='Length of PCR product, default: 150,2000.')

    parser.add_option('-p', '--proc',
                      dest='proc',
                      default="10",
                      type="int",
                      help='number of process.')

    parser.add_option('-m', '--seedmms',
                      dest='seedmms',
                      default="1",
                      type="int",
                      help='Mismatches in seed (can be 0 - 3, default: -n 1).')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='Prediction of off-targets PCR product.')

    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input_file is None:
        parser.print_help()
        print("Input file must be specified !!!")
        sys.exit(1)
    elif options.ref is None:
        parser.print_help()
        print("reference index (bowtie) must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("No output file provided !!!")
        sys.exit(1)
    return parser.parse_args()


TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")

degenerate_base = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

score_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
               "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}

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

adjust = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}

adjust_initiation = {"A": 2.8, "T": 2.8, "C": 1.82, "G": 1.82}

adjust_terminal_TA = 0.4

symmetry_correction = 0.4

base2bit = {"A": 0, "C": 1, "G": 2, "T": 3}

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


def RC(seq):
    return seq.translate(TRANS)[::-1]


def Penalty_points(length, GC, d1, d2):
    return log10((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)))

def nan_removing(pre_list):
    while np.nan in pre_list:
        pre_list.remove(np.nan)
    return pre_list
####################################################

def build_dict(Input):
    Input_dict = defaultdict(list)
    with open(Input, "r") as f:
        position_pattern_1 = re.compile('MD:Z:(\w+)')
        position_pattern = re.compile("[A-Z]?(\d+)")
        for i in f:
            i = i.strip().split("\t")
            primer = i[0]
            gene = i[2]
            primer_match_start = int(i[3]) - 1
            candidate_MD = nan_removing(i[11:])
            string = str('\t'.join(candidate_MD))
            position_1 = position_pattern_1.search(string).group(1)
            position = position_pattern.search(position_1[-2:]).group(1)
            if int(position) < 5:
                pass
            else:
                Input_dict[gene].append([primer_match_start, primer])
    return Input_dict


def closest(my_list, my_number1, my_number2):
    index_left = bisect_left(my_list, my_number1)
    # find the first element index in my_list which greater than my_number.
    if my_number2 > my_list[-1]:
        index_right = len(my_list) - 1  # This is index.
    else:
        index_right = bisect_left(my_list, my_number2) - 1
    return index_left, index_right


class off_targets(object):
    def __init__(self, primer_file, term_length=9, reference_file="",mismatch_num=1,
                 PCR_product_size="150,2000", outfile="", nproc=10):
        #  If an attribute in a Python class does not want to be accessed externally,
        #  we can start with a double underscore (__) when naming the attribute,
        #  Then the attribute cannot be accessed with the original variable name, making it private.
        #  If an attribute is marked with "__xxxx_" Is defined, then it can be accessed externally.
        self.nproc = nproc
        self.term_len = term_length
        self.primer_file = primer_file
        self.reference_file = reference_file
        self.outfile = outfile
        self.PCR_size = PCR_product_size
        self.resQ = Manager().Queue()
        self.mismatch_num = mismatch_num

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

    def get_term9(self):
        Output = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".term.fa")
        Input = self.primer_file
        length = self.term_len
        term_list = defaultdict(list)
        seq_ID = defaultdict(list)
        l = length
        with open(Input, "r") as f:
            for i in f:
                if i.startswith(">"):
                    value = i.strip()
                else:
                    key = i.strip()[-l:]
                    term_list[key].append(value)

        for k in term_list.keys():
            sequence = k
            Id = "|".join(term_list[k])
            expand_seq = self.degenerate_seq(sequence)
            if len(expand_seq) > 1:
                for j in range(len(expand_seq)):
                    ID = Id + "_" + str(j)
                    seq_ID[expand_seq[j]].append(ID)
            else:
                ID = Id + "_0"
                seq_ID[sequence].append(ID)
        with open(Output, "w") as fo:
            for seq in seq_ID.keys():
                fo.write('|'.join(seq_ID[seq]) + "\n" + seq + "\n")
        return seq_ID

    def bowtie_map(self):
        fa = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".term.fa")
        ref_index = self.reference_file
        out = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".sam")
        for_out = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".for.sam")
        rev_out = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".rev.sam")
        nproc = self.nproc
        if for_out.exists() and rev_out.exists():
            pass
        else:
            # os.system("bowtie2 -p 20 -f -N 0 -a -x {} -f -U {} -S {}".format(ref_index, fa, out))
            os.system(
                "bowtie -f -n {} -l 8 -a -p {} --best --strata {} {} -S {}".format(self.mismatch_num, nproc, ref_index,
                                                                                   fa, out))
            os.system("samtools view -@ {} -F 16 {} > {}".format(nproc, out, for_out))
            os.system("samtools view -@ {} -f 16 {} > {}".format(nproc, out, rev_out))

    def build_dict_run(self):
        sam_for_file = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".for.sam")
        sam_rev_file = Path(self.primer_file).parent.joinpath(Path(self.primer_file).stem).with_suffix(".rev.sam")
        pool = multiprocessing.Pool()
        # pool.apply_async(build_dict, kwds={sam_for_file, sam_rev_file})
        forward_dict, reverse_dict = pool.map_async(build_dict, (sam_for_file, sam_rev_file)).get()
        pool.close()
        pool.join()
        print("Number of genes with candidate primers: forward ==> {}; reverse ==> {}.".format(len(forward_dict),
                                                                                               len(reverse_dict)))
        target_gene = list(set(forward_dict.keys()).intersection(reverse_dict.keys()))
        print("Number of genes with candidate primer pairs: {}.".format(
            len(set(forward_dict.keys()).intersection(reverse_dict.keys()))))
        return target_gene, forward_dict, reverse_dict

    def PCR_prediction(self, gene, F_dict, R_dict):
        product_len = self.PCR_size.split(",")
        primer_F = dict(F_dict[gene])
        position_start = sorted(primer_F.keys())
        primer_R = dict(R_dict[gene])
        position_stop = sorted(primer_R.keys())
        if int(position_stop[0]) - int(position_start[-1]) > int(product_len[1]):
            pass
        elif int(position_stop[-1]) - int(position_start[0]) < int(product_len[0]):
            pass
        else:
            for start in range(len(position_start)):
                stop_index_start, stop_index_stop = closest(position_stop,
                                                                 position_start[start] + int(product_len[0]),
                                                                 position_start[start] + int(product_len[1]))
                if stop_index_start > stop_index_stop:  # caution: all(var) > stop_index_start in bisect_left,
                    # you need to stop immediately when stop_index_start > Product length
                    break
                else:
                    for stop in range(stop_index_start, stop_index_stop + 1):
                        distance = int(position_stop[stop]) - int(position_start[start]) + 1
                        if distance > int(product_len[1]):
                            break
                        elif int(product_len[0]) < distance < int(product_len[1]):
                            line = (gene, int(position_start[start]), int(position_stop[stop]),
                                    primer_F[position_start[start]], primer_R[position_stop[stop]], distance)
                            # off_target = {"Chrom (or Genes)": gene,
                            #               "Start": int(position_start[start]),
                            #               "Stop": int(position_stop[stop]),
                            #               "Primer_F": primer_F[position_start[start]],
                            #               "Primer_R": primer_R[position_stop[stop]],
                            #               "Product length": distance}
                            # out.append(off_target)
                            self.resQ.put(line)
                            # In multiple processes, each process has its own variable copy, so a variable in the
                            # main process is transferred to other processes for modification, and the result is
                            # still stored in that process. In fact, this variable in the main process is equivalent
                            # to not being modified. In order to synchronize the changes of other processes to the
                            # main process, you need to create variables that can be shared among multiple processes.
        self.resQ.put(None)

    def run(self):
        self.get_term9()
        self.bowtie_map()
        target_gene, forward_dict, reverse_dict = self.build_dict_run()
        p = ProcessPoolExecutor(self.nproc)
        for gene in target_gene:
            p.submit(self.PCR_prediction(gene, forward_dict, reverse_dict))
        # This will submit all tasks to one place without blocking, and then each
        # thread in the thread pool will fetch tasks.
        n = 0
        primer_pair_id = defaultdict(int)
        primer_id = defaultdict(int)
        #primer_reverse_id = defaultdict(int)
        with open(self.outfile, "w") as fo:
            headers = ["Chrom (or Genes)", "Start", "Stop", "Primer_F", "Primer_R", "Product length"]
            fo.write("\t".join(headers) + "\n")
            while n < len(target_gene):
                res = self.resQ.get()
                # The get method can read and delete an element from the queue. Similarly, the get method has two
                # optional parameters: blocked and timeout. If blocked is true (the default value) and timeout is
                # positive, no element is retrieved during the waiting time, and a Queue is thrown Empty exception.
                # If blocked is false, there are two cases. If a value of Queue is available, return the value
                # immediately. Otherwise, if the queue is empty, throw a Queue.Empty exception immediately.
                if res is None:
                    n += 1
                    continue
                primer_pair_id[res[3] + "\t" + res[4]] += 1
                primer_id[res[3]] += 1
                primer_id[res[4]] += 1
                fo.write("\t".join(map(str, res)) + "\n")
                # get results before shutdown. Synchronous call mode: call, wait for the return value, decouple,
                # but slow.
        p.shutdown()
        primer_pair_id_sort = sorted(primer_pair_id.items(), key=lambda x: x[1], reverse=True)
        primer_id_sort = sorted(primer_id.items(), key=lambda x: x[1], reverse=True)
        with open(self.outfile + ".pair.num", "w") as fo:
            fo.write("Primer_F\tPrimer_R\tPair_num\n")
            for k in primer_pair_id_sort:
                fo.write(k[0] + "\t" + str(k[1]) + "\n")
        with open(self.outfile + ".num", "w") as fo2:
            fo2.write("Primer\tnum\n")
            for k in primer_id_sort:
                fo2.write(k[0] + "\t" + str(k[1]) + "\n")

def main():
    options, args = argsParse()
    prediction = off_targets(primer_file=options.input_file, term_length=options.len, reference_file=options.ref,
                             PCR_product_size=options.size, mismatch_num=options.seedmms, outfile=options.out,
                             nproc=options.proc)
    prediction.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
