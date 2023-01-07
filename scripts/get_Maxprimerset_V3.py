#!/bin/python
"""
Get the final primerset.
construct primer_set by dimer_check, off_targets_prediction, deltaG_filter.
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
import multiprocessing
import re
import pickle
import os
import sys
import time
from bisect import bisect_left
from collections import defaultdict
from optparse import OptionParser
import re
import math
from operator import mul
from functools import reduce
import pandas as pd
import numpy as np
from pathlib import Path
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor
from itertools import product


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -o [output] \n \
                            Options: {-s [step] -m [T]}', version="%prog 0.0.4")
    parser.add_option('-c', '--core',
                      dest='core',
                      help='core primer file: core_primers.fa. Only fasta is accepted.')

    parser.add_option('-n', '--new',
                      dest='new',
                      help='new primer file: new_primers.fa. Only fasta is accepted.')

    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='reference file: bowtie index.')

    parser.add_option('-s', '--size',
                      dest='size',
                      default="150,2000",
                      help='Length of PCR product, default: [150,2000].')

    parser.add_option('-l', '--len',
                      dest='len',
                      default="9",
                      type="int",
                      help='Length of primer term, which is used for off-targets prediction.')

    parser.add_option('-p', '--proc',
                      dest='proc',
                      default="10",
                      type="int",
                      help='Number of process to launch.  default: 10.')

    parser.add_option('-m', '--model',
                      dest='model',
                      default="DO",
                      type="str",
                      help='D means finDimer; O means off-targets prediction. Default: DO.'
                           'Use D or O if you want fin[D]imer or predict [O]ff-targets.')

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
    return parser.parse_args()


# TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")
TRANS = str.maketrans("ATGC", "TACG")

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


################################################################
def RC(seq):
    return seq.translate(TRANS)[::-1]


def score_trans(sequence):
    return reduce(mul, [math.floor(score_table[x]) for x in list(sequence)])


def Penalty_points(length, GC, d1, d2):
    return math.log10((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)))


class Dimer(object):

    def __init__(self, core_primer_file, new_primer_set, outfile, nproc=10):
        self.nproc = nproc
        self.core_primer_file = core_primer_file
        self.new_primer_set = new_primer_set
        self.core_primer = self.parse_primers(core_primer_file)
        self.new_primer = self.parse_primers(new_primer_set)
        self.outfile = os.path.abspath(outfile)
        self.resQ = Manager().Queue()

    def parse_primers(self, Input):
        primer_dict = defaultdict(str)
        Output = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".primer_set")
        if Output.exists():
            primer_info = open(Output, "rb")
            primer_dict = pickle.load(primer_info)
        else:
            with open(Input, "r") as f:
                for i in f:
                    if i.startswith(">"):
                        name = i.strip()
                    else:
                        primer_dict[i.strip()] = name
            with open(Output, "wb") as fo:
                pickle.dump(primer_dict, fo)
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
        return ("".join(i) for i in product(*seq))

    def current_end(self, primer, adaptor="", num=5, length=14):
        primer_extend = adaptor + primer
        end_seq = []
        for i in range(num, (num + length)):
            s = primer_extend[-i:]
            if s:
                end_seq.extend(self.degenerate_seq(s))
        return end_seq

    #
    # def end_set(self, Input, num=5, length=14):
    #     Output = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".end")
    #     if Output.exists():
    #         end_info = open(Output, "rb")
    #         end_seq = pickle.load(end_info)
    #     else:
    #         end_seq = []
    #         with open(Input,"r") as fi:
    #             for line in fi:
    #                 if line.startswith("#"):
    #                     pass
    #                 elif line.startswith(">"):
    #                     pass
    #                 else:
    #                     for i in range(num, (num + length)):
    #                         s = line.strip()[-i:]
    #                         if s:
    #                             end_seq.extend(self.degenerate_seq(s))
    #         with open(Output, "wb") as fo:
    #             pickle.dump(end_seq, fo)
    #     return end_seq

    def deltaG(self, sequence):
        Delta_G_list = []
        for seq in self.degenerate_seq(sequence):
            Delta_G = 0
            for n in range(len(seq) - 1):
                i, j = base2bit[seq[n + 1]], base2bit[seq[n]]
                Delta_G += freedom_of_H_37_table[i][j] * H_bonds_number[i][j] + penalty_of_H_37_table[i][j]
            term5 = sequence[-2:]
            if term5 == "TA":
                Delta_G += adjust_initiation[seq[0]] + adjust_terminal_TA + symmetry_correction
            else:
                Delta_G += adjust_initiation[seq[0]] + symmetry_correction
            Delta_G_list.append(Delta_G)
        return round(max(Delta_G_list), 2)

    def dimer_check(self, primer, group):
        current_end_set = self.current_end(primer)
        current_end_list = sorted(list(current_end_set), key=lambda i: len(i), reverse=True)
        dimer = False
        if group == "core":
            primers = self.new_primer
        else:
            primers = self.core_primer
            # findimer between new_primer
            for ps in self.new_primer:
                for end in current_end_list:
                    for p in self.degenerate_seq(ps):
                        idx = p.find(RC(end))
                        if idx >= 0:
                            end_length = len(end)
                            end_GC = end.count("G") + end.count("C")
                            end_d1 = 0
                            end_d2 = len(p) - len(end) - idx
                            Loss = Penalty_points(
                                end_length, end_GC, end_d1, end_d2)
                            delta_G = self.deltaG(end)
                            if Loss > 3 or (delta_G < -5 and Loss >= 2):
                                line = (self.new_primer[primer], primer,
                                        end, delta_G, end_length, end_d1,
                                        end_GC, self.new_primer[ps],
                                        ps, end_d2, Loss
                                        )
                                self.resQ.put(line)
                                dimer = True
                                if dimer:
                                    break
                    if dimer:
                        dimer = False
                        break
            self.resQ.put(None)
        for ps in primers:
            # cross validation (findimer between new and core)
            for end in current_end_list:
                for p in self.degenerate_seq(ps):
                    idx = p.find(RC(end))
                    if idx >= 0:
                        end_length = len(end)
                        end_GC = end.count("G") + end.count("C")
                        end_d1 = 0
                        end_d2 = len(p) - len(end) - idx
                        Loss = Penalty_points(
                            end_length, end_GC, end_d1, end_d2)
                        delta_G = self.deltaG(end)
                        if Loss > 3 or delta_G < -5:
                            line = (primers[primer], primer,
                                    end, delta_G, end_length, end_d1,
                                    end_GC, primers[ps],
                                    ps, end_d2, Loss
                                    )
                            self.resQ.put(line)
                            dimer = True
                            if dimer:
                                break
                if dimer:
                    dimer = False
                    break
        self.resQ.put(None)  # Found a None, you can rest now

    def run(self):
        p = ProcessPoolExecutor(self.nproc)
        for primer in self.core_primer.keys():
            p.submit(self.dimer_check, primer=primer, group="core")
        for primer in self.new_primer.keys():
            p.submit(self.dimer_check, primer=primer, group="new")
        n = 0
        primer_id_sum = defaultdict(int)
        dimer_primer_id_sum = defaultdict(int)
        with open(self.outfile, "w") as fo:
            headers = ["Primer_ID", "Primer seq", "Primer end", "Delta G", "Primer end length", "End (distance 1)",
                       "End (GC)", "Dimer-primer_ID", "Dimer-primer seq", "End (distance 2)", "Loss"]
            fo.write("\t".join(headers) + "\n")
            while n < len(self.new_primer.keys()) + len(self.core_primer.keys()):
                res = self.resQ.get()
                if res is None:
                    n += 1
                    continue
                primer_id_sum[res[0]] += 1
                dimer_primer_id_sum[res[7]] += 1
                fo.write("\t".join(map(str, res)) + "\n")
        p.shutdown()
        with open(self.outfile + ".dimer_num", "w") as fo:
            fo.write("SeqName\tPrimer_ID\tDimer-primer_ID\tRowSum\n")
            for k in primer_id_sum.keys():
                p_id = primer_id_sum[k]
                d_id = dimer_primer_id_sum[k]
                RowSum = p_id + d_id
                fo.write("\t".join(map(str, [k, p_id, d_id, RowSum])) + "\n")


####################################################

class off_targets(object):
    def __init__(self, core_primer_file, new_primer_set, term_length=9, reference_file="",
                 PCR_product_size="150,2000", outfile="", nproc=10):
        self.nproc = nproc
        self.term_len = term_length
        self.core_primer_file = core_primer_file
        self.new_primer_set = new_primer_set
        self.reference_file = reference_file
        self.outfile = outfile
        self.PCR_size = PCR_product_size
        self.resQ = Manager().Queue()

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

    def get_term(self, Input):
        Output = Path(Input).parent.joinpath(Path(Input.stem)).with_suffix(".term.fa")
        if Output.exists():
            pass
        else:
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

    def bowtie_map(self, Input):
        fa = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".term.fa")
        ref_index = self.reference_file
        out = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".sam")
        for_out = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".for.sam")
        rev_out = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".rev.sam")
        nproc = self.nproc
        if for_out.exists() and rev_out.exists():
            pass
        else:
            # os.system("bowtie2 -p 20 -f -N 0 -a -x {} -f -U {} -S {}".format(ref_index, fa, out))
            os.system(
                "bowtie -f -n 1 -l 9 -a -p {} --best --strata {} {} -S {}".format(nproc, ref_index, fa, out))
            os.system("samtools view -@ {} -F 16 {} > {}".format(nproc, out, for_out))
            os.system("samtools view -@ {} -f 16 {} > {}".format(nproc, out, rev_out))

    def build_dict(self, Input):
        Output = Path(Input).parent.joinpath(Input.stem).with_suffix(".gene_position")
        if Output.exists():
            with open(Output, "rb") as fi:
                Input_dict = pickle.load(fi)
        else:
            Input_dict = defaultdict(list)
            with open(Input, "r") as f:
                for i in f:
                    i = i.strip().split("\t")
                    primer = i[0]
                    gene = i[2]
                    primer_match_start = int(i[3]) - 1
                    Input_dict[gene].append([primer_match_start, primer])
            with open(Output, "wb") as fo:
                pickle.dump(Input_dict, fo)
        return Input_dict

    @staticmethod
    def closest(my_list, my_number1, my_number2):
        index_left = bisect_left(my_list, my_number1)
        # find the first element index in my_list which greater than my_number.
        if my_number2 > my_list[-1]:
            index_right = len(my_list) - 1  # This is index.
        else:
            index_right = bisect_left(my_list, my_number2) - 1
        return index_left, index_right

    def build_dict_run(self, Input):
        sam_for_file = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".for.sam")
        sam_rev_file = Path(Input).parent.joinpath(Path(Input).stem).with_suffix(".rev.sam")
        pool = multiprocessing.Pool()
        forward_dict, reverse_dict = pool.map_async(self.build_dict, (sam_for_file, sam_rev_file)).get()
        pool.close()
        pool.join()
        print("Number of genes with candidate primers ({}): forward ==> {}; reverse ==> {}."
              .format(Input, len(forward_dict), len(reverse_dict)))
        return [forward_dict, reverse_dict]

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
                stop_index_start, stop_index_stop = self.closest(position_stop,
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
                            self.resQ.put(line)
        self.resQ.put(None)

    def run(self):
        pool_term = multiprocessing.Pool()
        pool_term.map_async(self.get_term, (self.core_primer_file, self.new_primer_set)).get()
        pool_term.close()
        pool_term.join()

        pool_dict = multiprocessing.Pool()
        core_list, new_list = pool_dict.map_async(self.build_dict_run,
                                                  (self.core_primer_file, self.new_primer_set)).get()
        pool_dict.close()
        pool_dict.join()

        core_forward_dict, core_reverse_dict = core_list[0], core_list[1]
        new_forward_dict, new_reverse_dict = new_list[0], new_list[1]

        target_gene_coreF_new_R = list(set(core_forward_dict.keys()).intersection(new_reverse_dict.keys()))
        target_gene_coreR_new_F = list(set(core_reverse_dict.keys()).intersection(new_forward_dict.keys()))
        target_gene_newF_new_R = list(set(new_forward_dict.keys()).intersection(new_reverse_dict.keys()))
        print("Number of genes with candidate primer pairs of coreF_new_R: {}.".format(len(target_gene_coreF_new_R)))
        print("Number of genes with candidate primer pairs of coreR_new_F: {}.".format(len(target_gene_coreR_new_F)))
        print("Number of genes with candidate primer pairs of newF_new_R: {}.".format(len(target_gene_newF_new_R)))
        p = ProcessPoolExecutor(self.nproc)
        primer_forward_id = defaultdict(int)
        primer_reverse_id = defaultdict(int)
        for gene1 in target_gene_coreF_new_R:
            p.submit(self.PCR_prediction(gene1, core_forward_dict, new_reverse_dict))
        for gene2 in target_gene_coreR_new_F:
            p.submit(self.PCR_prediction(gene2, new_forward_dict, core_reverse_dict))
        for gene3 in target_gene_newF_new_R:
            p.submit(self.PCR_prediction(gene3, new_forward_dict, new_reverse_dict))
        n = 0
        with open(self.outfile, "w") as fo:
            headers = ["Chrom (or Genes)", "Start", "Stop", "Primer_F", "Primer_R", "Product length"]
            fo.write("\t".join(headers) + "\n")
            while n < len(target_gene_coreF_new_R) + len(target_gene_coreR_new_F) + len(target_gene_newF_new_R):
                res = self.resQ.get()
                if res is None:
                    n += 1
                    continue
                primer_forward_id[res[3]] += 1
                primer_reverse_id[res[4]] += 1
                fo.write("\t".join(map(str, res)) + "\n")
        p.shutdown()
        with open(self.outfile + ".num", "w") as fo:
            fo.write("SeqName\tPrimer_ID\tDimer-primer_ID\tRowSum\n")
            for k in primer_forward_id.keys():
                p_id = primer_forward_id[k]
                d_id = primer_reverse_id[k]
                RowSum = p_id + d_id
                fo.write("\t".join(map(str, [k, p_id, d_id, RowSum])) + "\n")


def Dimer_main():
    # options, args = argsParse()
    Dimer_out = options.out + ".dimer"
    dimer_app = Dimer(core_primer_file=options.core, new_primer_set=options.new, outfile=Dimer_out,
                      nproc=options.proc)
    dimer_app.run()


def Offtargets_main():
    # options, args = argsParse()
    off_targets_out = options.out + ".offtargets"
    prediction = off_targets(core_primer_file=options.core, new_primer_set=options.new, term_length=options.len,
                             reference_file=options.ref, PCR_product_size=options.size, outfile=off_targets_out,
                             nproc=options.proc)
    prediction.run()


if __name__ == "__main__":
    e1 = time.time()
    options, args = argsParse()
    if re.search("D", options.model):
        Dimer_main()
    elif re.search("O", options.model):
        Offtargets_main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
