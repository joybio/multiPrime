#!/bin/python
""" Stastic of primer pair coverage by multiPrime prediction (top N) and BWT (total) validation."""
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
    parser = OptionParser('Usage: Stastics of primer pair coverage by multiPrime prediction (top N) and BWT.'
			'%prog -i [input] -r [Cluster_0_20723.fa] -c [.pair.num] -o [output]')

    parser.add_option('-i', '--input',
                      dest='input',
                      help='input file: candidate primer.txt. Results of get_multiPrime.py, which located in the dir: Cluster_cprimer.')

    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='reference file, all fasta clusters, for example: Cluster_0_20723.fa.')

    parser.add_option('-c', '--coverage',
                      dest='coverage',
                      help='BWT coverage results.')

    parser.add_option('-s', '--step',
                      dest='step',
                      default="5",
                      type="int",
                      help='Step length. Step between Primer_1_F : Primer_2_F. Default: 5')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='Output files.')

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

def position_coverage(Input, step):
    primer_coverage_dict = {}
    with open(Input, "r") as f:
            for i in f:
                if i.startswith("#"):
                    pass
                else:
                    i = i.strip().split("\t")
                    # primer = i[0].split("/")[-1].strip(".candidate.primers.txt")
                    n = 0
                    while n < len(i)-3:
                            coverage = i[n + 3].split(":")[2]
                            position = i[n + 5]
                            primer_coverage_dict[position] = coverage
                            n += step
    return primer_coverage_dict

def get_number(Input):
    from itertools import (takewhile, repeat)
    buffer = 1024 * 1024
    with open(Input, encoding="utf-8") as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        seq_number = int(sum(buf.count("\n") for buf in buf_gen) / 2)
    return seq_number

def correlation(Input, dict, seq_number, out):
    output = open(out, "w")
    header = ["Primer_F", "Primer_R", "coverage_of_multiPrime", "estimate_coverage_by BWT"]
    output.write("\t".join(header) + "\n")
    with open(Input) as f:
        for i in f:
            if i.startswith("Primer_F"):
                pass
            else:
                i = i.strip().split("\t")
                start = i[0].split("_")[-2]
                stop = i[1].split("_")[-2]
                key = str(start) + ":" + str(stop)
                if key in dict.keys():
                    primer_info = [i[0], i[1], dict[key], round(int(i[3])/seq_number,4)]
                    output.write("\t".join(map(str, primer_info)) + "\n")
    f.close()
    output.close()



def main():
    (options, args) = argsParse()
    # build dict
    primer_txt = options.input
    step = options.step
    primer_txt_dict = position_coverage(primer_txt, step)
    # number of input sequences
    BWT_results = options.coverage
    seq_number = get_number(options.ref)
    print(seq_number)
    # output
    correlation(BWT_results, primer_txt_dict, seq_number, options.out)


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))



