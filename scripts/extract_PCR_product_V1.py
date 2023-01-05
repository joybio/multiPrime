#!/bin/python
"""extract perfect PCR product by primers (degenerate is also ok!) with raw input fasta file."""

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
import time
from collections import defaultdict
from itertools import product
from multiprocessing import Manager
from optparse import OptionParser
from pathlib import Path
import os
from operator import mul
from functools import reduce
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from itertools import (takewhile, repeat)


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -p [primerF,primerR] -f [format] -o [output]', version="%prog 0.0.2")
    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='referebce file: template fasta or reference fasta.')

    parser.add_option('-i', '--input',
                      dest='input',
                      help='Primer file. One of the followed three types:\n '
                           'final_maxprimers_set.xls \n primer.fa \n primer_F,primer_R.')

    parser.add_option('-f', '--format',
                      dest='format',
                      help='Format of primer file: xls or fa or seq; default: xls. \n '
                           'xls: final_primer_set.xls, output of multiPrime. \n'
                           'fa: fasta format. \n'
                           'seq: sequence format, comma seperate. e.g. primer_F,Primer_R.')

    parser.add_option('-o', '--out',
                      dest='out',
                      default="PCR_product",
                      help='Output_dir. default: PCR_product.')

    parser.add_option('-p', '--process',
                      dest='process',
                      default="10",
                      type="int",
                      help='Number of process to launch.  default: 10.')

    parser.add_option('-s', '--stast',
                      dest='stast',
                      default="Coverage.xls",
                      help='Stast information: number of coverage and total. default: Coverage.xls')
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.ref is None:
        parser.print_help()
        print("Input (reference) file must be specified !!!")
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Primer file or sequence must be specified !!!")
        sys.exit(1)
    elif options.format is None:
        parser.print_help()
        print("Primer file format must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("No output file provided !!!")
        sys.exit(1)
    return parser.parse_args()


degenerate_base = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

score_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
               "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}

TRANS = str.maketrans("ATGC", "TACG")


def RC(seq):
    return seq.translate(TRANS)[::-1]


class Product(object):
    def __init__(self, primer_file="", output_file="", ref_file="", file_format="fa", coverage="", nproc=10):
        self.nproc = nproc
        self.primers_file = primer_file
        self.ref_file = ref_file
        self.output_file = Path(output_file)
        self.file_format = file_format
        self.primers = self.parse_primers()
        self.coverage = coverage
        self.resQ = Manager().Queue()

    def md_out_File(self):
        if self.output_file.exists():
            pass
        else:
            os.system("mkdir -p {}".format(self.output_file))

    def parse_primers(self):
        res = {}
        if self.file_format == "seq":
            primers = self.primers_file.split(",")
            res["PCR_info"] = [primers[0], primers[1]]
        else:
            with open(self.primers_file, "r") as f:
                if self.file_format == "xls":
                    for i in f:
                        if i.startswith("#"):
                            pass
                        else:
                            i = i.strip().split("\t")
                            cluster_id = i[0].split("/")[-1].split(".")[0]
                            primer_F = i[2]
                            primer_R = i[3]
                            start = i[6].split(":")[0]
                            stop = i[6].split(":")[1]
                            key = cluster_id + "_" + str(start) + "_F_" + cluster_id + "_" + str(stop)
                            res[key] = [primer_F, primer_R]
                elif self.file_format == "fa":
                    primer_info = pd.read_table(f, header=None)
                    for idx, row in primer_info.iterrows():
                        if idx % 4 == 0:
                            primer_F_info = row[0].lstrip(">")
                        elif idx % 4 == 1:
                            primer_F = row[0]
                        elif idx % 4 == 2:
                            key = primer_F_info + "_" + row[0].lstrip(">")
                        elif idx % 4 == 3:
                            primer_R = row[0]
                            res[key] = [primer_F, primer_R]
        return res

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
        # return ("".join(i) for i in product(*seq)) # This is a generator, just once when iteration
        # d = [x for x in range(12)]
        # g = (x for i in range(12))
        # The result of list derivation returns a list, and the tuple derivation returns a generator
        return ["".join(i) for i in product(*seq)]

    def get_PCR_PRODUCT(self, primerinfo, F, R, ref):
        Fseq = self.degenerate_seq(F)
        Rseq = self.degenerate_seq(R)
        product_dict = {}
        Non_targets_dict = {}
        with open(ref, "r") as r:
            for i in r:
                if i.startswith(">"):
                    key = i.strip()
                else:
                    value = ''
                    for sequence in Fseq:
                        if re.search(sequence, i):
                            line = i.split(sequence)
                            Product = sequence + line[1]
                            for sequence2 in Rseq:
                                if re.search(RC(sequence2), Product):
                                    Product = Product.split(RC(sequence2))
                                    value = Product[0].strip() + RC(sequence2)
                                    break
                            if value:
                                break
                    if value:
                        product_dict[key] = value
                    else:
                        Non_targets_dict[key] = i.strip()
        self.resQ.put([primerinfo, product_dict, Non_targets_dict, F, R])
        self.resQ.put(None)

    def run(self):
        self.md_out_File()
        proc = ProcessPoolExecutor(self.nproc)
        for primer in self.primers.keys():
            proc.submit(self.get_PCR_PRODUCT, primer, self.primers[primer][0], self.primers[primer][1], self.ref_file)
            #  This will submit all tasks to one place without blocking, and then each
            #  thread in the thread pool will fetch tasks
        n = 0
        Product_seq_id = set()
        non_Product_seq_id = set()
        while n < len(self.primers):
            res = self.resQ.get()
            if res is None:
                n += 1
                continue
            PCR_product = Path(self.output_file).joinpath(res[0]).with_suffix(".PCR.product.fa")
            PCR_non_product = Path(self.output_file).joinpath(res[0]).with_suffix(
                ".non_PCR.product.fa")
            with open(self.coverage, "a+") as c:
                c.write(
                    "Number of Product/non_Product, primer-F and primer-R: {}"
                    "\t{}\t{}\t{}\t{}\n".format(res[0], len(res[1].keys()), len(res[2].keys()), res[3], res[4]))
                with open(PCR_product, "w") as p:
                    for result in res[1].keys():
                        Product_seq_id.add(result)
                        p.write(result + "\n" + res[1][result] + "\n")
                with open(PCR_non_product, "w") as np:
                    for result2 in res[2].keys():
                        non_Product_seq_id.add(result2)
                        np.write(result2 + "\n" + res[2][result2] + "\n")
        proc.shutdown()
        # After I run the main, I don't care whether the sub thread is alive or dead. With this parameter, after all
        # the sub threads are executed, the main function is executed.
        # get results after shutdown. Asynchronous call mode: only call, unequal return values, coupling may exist,
        # but the speed is fast.
        buffer = 1024 * 1024
        with open(self.ref_file, encoding="utf-8") as f:
            buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
            seq_number = int(sum(buf.count("\n") for buf in buf_gen) / 2)
        with open(self.coverage, "a+") as c:
            c.write(
                "Total number of sequences:\t{}\n"
                "Coveraged number of sequence:\t{}\n"
                "Rate of coverage:\t> {}\n".format(seq_number, len(Product_seq_id),
                                                   round(float(len(Product_seq_id)) / seq_number, 2)))
        c.close()


def main():
    (options, args) = argsParse()
    results = Product(primer_file=options.input, output_file=options.out, ref_file=options.ref,
                      file_format=options.format, coverage=options.stast, nproc=options.process)
    results.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
