#!/bin/python

__date__ = "2022-10-8"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

import json
import os

import math
import time
from functools import reduce
from math import log10
from itertools import product
from multiprocessing import Manager
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor
from operator import mul
from statistics import mean
from optparse import OptionParser
import sys
import numpy as np
from itertools import repeat

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

import time
from multiprocessing import Manager
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from optparse import OptionParser
import sys


def argsParse():
    parser = OptionParser('Usage: %prog -i input -o output -p 10')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file.')

    parser.add_option('-p', '--nproc',
                      dest='nproc',
                      default="20",
                      type="int",
                      help="Number of process to launch. Default: 20.")

    parser.add_option('-t', '--threshold',
                      dest='threshold',
                      default="0",
                      type="int",
                      help="Those cluster sequence number less than threshold will try to merge to other clasters by ANI."
                           "Default: 0. "
                           "set -t 0, then all cluster will process."
                           "set -t 1, then no cluster will process.")

    parser.add_option('-o', '--out',
                      dest='out',
                      default="history.txt",
                      type='str',
                      help='Output file: history.txt.')
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


class merge_clstr(object):

    def __init__(self, inputfile="", output="", threshold=20, nproc=10):
        self.cluster_file = inputfile
        self.threshold = threshold
        self.work_dir = self.parse_work_dir()
        self.out = output
        self.cluster = self.parse_cluster()
        self.nproc = nproc
        self.resL = Manager().list()

    def parse_work_dir(self):
        return self.cluster_file.rstrip("cluster.txt") + "Clusters_fa"

    def parse_cluster(self):
        cluster_dict = defaultdict(int)
        with open(self.cluster_file, "r") as f:
            for i in f:
                if i.startswith("#"):
                    pass
                else:
                    i = i.strip().split("\t")
                    cluster_dict[i[0]] = int(i[1])
        cluster_dict_sort = sorted(cluster_dict.items(), key=lambda x: x[1], reverse=True)
        return cluster_dict_sort

    def merge_cluster_app(self, position, cluster_info):
        processing_ID = self.work_dir + "/" + cluster_info[position][0] + "_" + str(cluster_info[position][1])
        processing_file = processing_ID + ".txt"
        for i in range(len(self.cluster)):
            ref_ID = self.work_dir + "/" + self.cluster[i][0] + "_" + str(self.cluster[i][1])
            ref_file = ref_ID + ".txt"
            out_file = processing_ID + "_" + self.cluster[i][0] + "_" + str(self.cluster[i][1])
            if cluster_info[position][1] < self.cluster[i][1]:
                os.system("fastANI --ql {} --rl {} -o {}".format(processing_file, ref_file, out_file))
                ###### check primer #####
                size = os.path.getsize(out_file)
                if (size != 0):
                    print("Ref: {}, merge: {}".format(self.cluster[i], cluster_info[position]))
                    if self.cluster[i][1] > cluster_info[position][1]:
                        merge_info = [ref_ID, processing_ID]
                    else:
                        merge_info = [processing_ID, ref_ID]
                    self.resL.append(merge_info)
                    os.remove(out_file)
                    break
                else:
                    os.remove(out_file)
                    pass

    def run(self):
        p = ProcessPoolExecutor(self.nproc)
        cluster_info = self.cluster[::-1].copy()
        for position in range(len(self.cluster)):
            if self.threshold == 0:
                p.submit(self.merge_cluster_app, position, cluster_info)
            elif self.threshold == 1:
                break
            else:
                print("Abundunce: {}".format(cluster_info[position][1]))
                if cluster_info[position][1] <= self.threshold:
                    print("Position: {}".format(position))
                    p.submit(self.merge_cluster_app, position, cluster_info)
        p.shutdown()
        merge_dict = defaultdict(list)
        with open(self.out, "w") as f:
            for i in self.resL:
                f.write("\t".join(map(str, i)) + "\n")
                merge_dict[i[0]].append(i[1])
        f.close()
        for i in self.cluster:
            ANI_seq_ID = self.work_dir + "/" + i[0] + "_" + str(i[1])
            os.system("rm -rf {}".format(ANI_seq_ID))

        for ref in merge_dict.keys():
            ref_cluster_info = ref.split("_")
            ref_cluster = '_'.join(ref_cluster_info[:-1])
            ref_cluster_number = int(ref_cluster_info[-1])
            ref_fa = ref + ".fa"
            ref_tfa = ref + ".tfa"
            ref_txt = ref + ".txt"
            for sub in merge_dict[ref]:
                sub_cluster_info = sub.split("_")
                sub_cluster_number = sub_cluster_info[-1]
                ref_cluster_number += int(sub_cluster_number)
                sub_fa = sub + ".fa"
                sub_tfa = sub + ".tfa"
                sub_txt = sub + ".txt"
                os. system("cat {} >> {}".format(sub_fa, ref_fa))
                os. system("cat {} >> {}".format(sub_tfa, ref_tfa))
                os. system("cat {} >> {}".format(sub_txt, ref_txt))
                os.remove(sub_fa)
                os.remove(sub_tfa)
                os.remove(sub_txt)
            final_cluster = ref_cluster + "_" + str(ref_cluster_number)
            final_cluster_fa = final_cluster + ".fa"
            final_cluster_tfa = final_cluster + ".tfa"
            final_cluster_txt = final_cluster + ".txt"
            # os.rename(ref, final_cluster)
            os.rename(ref_fa, final_cluster_fa)
            os.rename(ref_tfa, final_cluster_tfa)
            os.rename(ref_txt, final_cluster_txt)

def main():
    options, args = argsParse()
    merge_clstr_app = merge_clstr(inputfile=options.input, output=options.out, threshold=options.threshold)
    merge_clstr_app.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
