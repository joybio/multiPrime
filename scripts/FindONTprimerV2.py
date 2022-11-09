#!/bin/python

__date__ = "2022-10-8"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

import difflib
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from multiprocessing import Manager
from optparse import OptionParser
import subprocess as sp
import gzip

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


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -v ["IVC"] -s [primer set] -p [20] -o  [output].')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: degeprimer out.')

    parser.add_option('-s', '--set',
                      dest='set',
                      help='primer set file.')

    parser.add_option('-p', '--nproc',
                      dest='nproc',
                      default=10,
                      type="int",
                      help='Primer set file. option. Default: 10')

    parser.add_option('-l', '--len',
                      dest='len',
                      default=18,
                      type="int",
                      help='Primer length. Default: 18')

    parser.add_option('-m', '--min_ident',
                      dest='min_ident',
                      default=0.6,
                      type="float",
                      help='min identity. Default: 0.6')

    parser.add_option('-f', '--format',
                      dest='format',
                      default="fq",
                      type="str",
                      help='Input format, fq or fa, Default: fq')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='Output file: candidate primers. e.g. [*].candidate.primers.txt.')
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


degenerate_base = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}


# def build_blast_database_if_needed(seqs):  # 判断blastn所需要的db文件是否存在
#     if not os.path.exists(seqs + '.nin'):
#         with open(os.devnull, 'w') as devnull:
#             sp.run('makeblastdb -dbtype nucl -in ' + seqs, stdout=devnull, shell=True, check=True)
#

# def iter_count(file_name):
#     from itertools import (takewhile, repeat)
#     buffer = 1024 * 1024
#     if file_name.endswith("gz"):
#         with gzip.open(file_name, "rb") as f:
#             buf_gen = takewhile(lambda x: x, (f.read(buffer).decode() for _ in repeat(None)))
#             # 读取到缓冲区
#             return int(sum(buf.count("\n") for buf in buf_gen) / 4)
#     else:
#         with open(file_name, "r") as f:
#             buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
#             return int(sum(buf.count("\n") for buf in buf_gen) / 4)


class ONTprimer(object):
    def __init__(self, input_file, primer_file, outfile="",  primer_len=18, min_ident=0.8, Input_format="fq", nproc=10):
        self.input_file = input_file
        self.primer_file = primer_file
        self.expand_dict = self.get_expand()
        self.Input_format = Input_format
        self.primer_len = primer_len
        self.outfile = outfile
        self.min_ident = min_ident
        self.nproc = nproc
        self.resQ = Manager().Queue()

    # def run_blastn(self, query):
    #     db = self.primer_file.strip("fa") + "expand.fa"
    #     min_identity = self.min_ident
    #     build_blast_database_if_needed(db)
    #
    #     cmd = 'blastn -db {} -query {}'.format(db, query)  # .format()感觉比用%s，%d更好
    #     cmd += " -outfmt 6 -evalue 1e-5 -max_target_seqs 1"
    #     cmd += ' -perc_identity {}'.format(min_identity)
    #
    #     f = os.popen(cmd)
    #     blast_hits = f.read().split("\n")
    #     f.close()
    #     return blast_hits

    @staticmethod
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
        return ["".join(i) for i in product(*seq)]

    def get_expand(self):
        expand_dict = {}
        expand_file = self.primer_file.strip("fa") + "expand.fa"
        with open(expand_file, "w") as f:
            with open(self.primer_file, "r") as p:
                for i in p:
                    if i.startswith(">"):
                        ID = i.strip()
                    else:
                        seq = self.degenerate_seq(i.strip())
                        for j in range(len(seq)):
                            f.write(ID + " | " + str(j) + "\n" + seq[j] + "\n")
                            expand_dict[seq[j]] = ID + " | " + str(j)
        p.close()
        f.close()
        return expand_dict

    def get_primer(self, seq_list):
        expand_list = list(self.expand_dict.keys())
        primer_F_R_list = []
        for p in range(len(seq_list)):
            expand_ratio_list = []
            seq = seq_list[p]
            for j in range(len(expand_list)):
                ratio = round(difflib.SequenceMatcher(None, seq, expand_list[j]).ratio(), 2)
                expand_ratio_list.append(ratio)
            max_ratio = max(expand_ratio_list)
            max_ratio_idx = expand_ratio_list.index(max_ratio)
            if max_ratio > self.min_ident:
                primer = self.expand_dict[expand_list[max_ratio_idx]].split(" | ")[0]
                # 字符"?"重复前面一个匹配字符零次或者一次.
            else:
                primer = "NA"
            primer_F_R_list.append(primer)
        self.resQ.put("\t".join(primer_F_R_list))
        # 给get_primer传输一个空的数据, 使get_primer模块运行结束.
        # Found a None, you can rest now

    def run(self):
        p = ProcessPoolExecutor(self.nproc)
        seq_number = 0
        if self.Input_format == "fq":
            if self.input_file.endswith("gz"):
                f = gzip.open(self.input_file, 'rb')
                for idx, line in enumerate(f.readlines()):  # 按行进行读取
                    if idx % 4 == 1:
                        s = line.decode().strip()  # 读取之后要进行解码
                        # print(s)
                        Primer_F_R = [s[:self.primer_len], s[-self.primer_len:]]
                        # print(Primer_F_R)
                        p.submit(self.get_primer, Primer_F_R)
                        seq_number += 1
                        # print(seq_number)
                        #  This will submit all tasks to one place without blocking, and then each
                        #  thread in the thread pool will fetch tasks
            else:
                f = open(self.input_file, 'r')
                for idx, line in enumerate(f.readlines()):  # 按行进行读取
                    if idx % 4 == 1:
                        Primer_F_R = [line[:self.primer_len], line[-self.primer_len:]]
                        p.submit(self.get_primer, Primer_F_R)
                        seq_number += 1
                        # print(seq_number)
        elif self.Input_format == "fa":
            if self.input_file.endswith("gz"):
                f = gzip.open(self.input_file, 'rb')
                for idx, line in enumerate(f.readlines()):  # 按行进行读取
                    if idx % 2 == 1:
                        s = line.decode().strip()  # 读取之后要进行解码
                        # print(s)
                        Primer_F_R = [s[:self.primer_len], s[-self.primer_len:]]
                        # print(Primer_F_R)
                        p.submit(self.get_primer, Primer_F_R)
                        seq_number += 1
                        # print(seq_number)
                        #  This will submit all tasks to one place without blocking, and then each
                        #  thread in the thread pool will fetch tasks
            else:
                f = open(self.input_file, 'r')
                for idx, line in enumerate(f.readlines()):  # 按行进行读取
                    if idx % 2 == 1:
                        Primer_F_R = [line[:self.primer_len], line[-self.primer_len:]]
                        p.submit(self.get_primer, Primer_F_R)
                        seq_number += 1
                        # print(seq_number)
        else:
            print("Please command line! only fa or fq is accepted!")
            sys.exit(1)
        n = 0
        primer_pair_sum = defaultdict(int)
        # seq_number = iter_count(self.input_file)
        print(seq_number)
        while n < seq_number:
            res = self.resQ.get()
            # The get method can read and delete an element from the queue. Similarly, the get method has two
            # optional parameters: blocked and timeout. If blocked is true (the default value) and timeout is
            # positive, no element is retrieved during the waiting time, and a Queue is thrown Empty exception.
            # If blocked is false, there are two cases. If a value of Queue is available, return the value
            # immediately. Otherwise, if the queue is empty, throw a Queue.Empty exception immediately.
            # if res is None:
            #     n += 1
            #     continue
            # # print(res)
            primer_pair_sum[res] += 1
            n += 1
            # print(n)
        # try:
        #     res = self.resQ.get()
        # except self.resQ.empty():
        #     break
        print(len(primer_pair_sum))
        count_sort = sorted(primer_pair_sum.items(), key=lambda x: x[1], reverse=True)
        with open(self.outfile + ".num", "w") as fo:
            headers = ["Primer_F", "Primer_R", "Number"]
            fo.write("\t".join(headers) + "\n")
            for i in count_sort:
                fo.write("\t".join(map(str, i)) + "\n")
        fo.close()
        #  get results before shutdown. Synchronous call mode: call, wait for the return value, decouple, but slow.
        p.shutdown()


def main():
    options, args = argsParse()
    ONTprimer_app = ONTprimer(options.input, options.set, outfile=options.out, primer_len=options.len,
                              Input_format=options.format, min_ident=options.min_ident, nproc=options.nproc)
    ONTprimer_app.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
