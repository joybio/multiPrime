#!/bin/python
# coding:utf-8
"""
Output off-target PCR results. 
"""
__date__ = "2022-9-27"
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


import sys
from collections import defaultdict
import os
import multiprocessing
from multiprocessing import Process
import time
import pandas as pd
from optparse import OptionParser
from pathlib import Path
import math
from functools import reduce
from operator import mul  #
from bisect import bisect_left


# Path(path).parent, Path(path).name, Path(path).suffix, Path(path).stem, Path(path).iterdir(), Path(path).joinpath()
# Path(path).is_absolute(), Path(path).is_dir(), Path(path).is_file(), Path(path).exists(), Path(path).with_suffix


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -r [bowtie index] -l [150,2000] -p [10]-o [output]')

    parser.add_option('-i', '--input',
                      dest='input',
                      help='input file: primer.fa.')

    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='reference file: bowtie index.')

    parser.add_option('-l', '--len',
                      dest='len',
                      default="150,2000",
                      help='Length of PCR product, default: 150,2000.')

    parser.add_option('-p', '--proc',
                      dest='proc',
                      default="10",
                      type="int",
                      help='number of process.')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='Prediction of off-targets PCR product.')

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
        print("reference index (bowtie) must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("No output file provided !!!")
        sys.exit(1)
    return parser.parse_args()


degenerate_pair = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

degenerate_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
                    "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}


def score_trans(sequence):
    return reduce(mul, [math.floor(degenerate_table[x]) for x in list(sequence)])


def dege_trans(sequence):
    expand_seq = [sequence.upper()]
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


def get_term9(fa, out):
    # primer_list = []
    with open(fa, "r") as f:
        with open(out, "w") as o:
            for i in f:
                if i.startswith(">"):
                    # primer = i.lstrip(">").rstrip()
                    # primer_list.append(primer)
                    o.write(i)
                else:
                    i = i.strip()[-9:]
                    o.write(i + "\n")
    # return primer_list


def bowtie_map(fa, ref_index, out, for_out, rev_out, number):
    # os.system("bowtie2 -p 20 -f -N 0 -a -x {} -f -U {} -S {}".format(ref_index, fa, out))
    os.system(
        "/share/data3/yangjunbo/miniconda3/bin/bowtie -f -n 1 -l 9 -a -p {} --best --strata {} {} -S {}".format(
            number, ref_index, fa, out))
    os.system("samtools view -@ {} -F 16 {} > {}".format(number, out, for_out))
    os.system("samtools view -@ {} -f 16 {} > {}".format(number, out, rev_out))


def build_dict(Input):
    Input_dict = defaultdict(list)
    with open(Input, "r") as f:
        for i in f:
            i = i.strip().split("\t")
            primer = i[0]
            gene = i[2]
            primer_match_start = int(i[3]) - 1
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


def PCR_prediction(target_gene, F_dict, R_dict, out):
    for gene in target_gene:
        primer_F = dict(F_dict[gene])
        position_start = sorted(primer_F.keys())
        primer_R = dict(R_dict[gene])
        position_stop = sorted(primer_R.keys())
        # print(position_stop, position_start)
        if int(position_stop[0]) - int(position_start[-1]) > int(product_len[1]):
            # print("min product length is: {}; next ...".format(int(position_stop[0]) - int(position_start[-1])))
            pass
        elif int(position_stop[-1]) - int(position_start[0]) < int(product_len[0]):
            # print("max product length is: {}; next ...".format(int(position_stop[-1]) - int(position_start[0])))
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
                            off_target = {"Chrom (or Genes)": gene,
                                          "Start": int(position_start[start]),
                                          "Stop": int(position_stop[stop]),
                                          "Primer_F": primer_F[position_start[start]],
                                          "Primer_R": primer_R[position_stop[stop]],
                                          "Product length": distance}
                            # 就相当于先向服务器请求数据 (此处不需要请求), 再向服务器传送修改后的数据.
                            out.append(off_target)


def get_dict_key(d_dict, value):
    d_key = list(d_dict.keys())[list(d_dict.values()).index(value)]
    return d_key


def fa_dege_trans(fa, output):
    seq_ID = defaultdict(list)
    with open(fa, "r") as dege_fa:
        for i in dege_fa:
            if i.startswith(">"):
                Id = i.strip()
            else:
                sequence = i.strip()
                expand_seq = dege_trans(sequence)
                if len(expand_seq) > 1:
                    for j in range(len(expand_seq)):
                        ID = Id + "_" + str(j)
                        seq_ID[expand_seq[j]].append(ID)
                else:
                    ID = Id + "_0"
                    seq_ID[sequence].append(ID)
    primer_list = []
    with open(output, "w") as out:
        for seq in seq_ID.keys():
            primer_list.append('|'.join(seq_ID[seq]))
            out.write('|'.join(seq_ID[seq]) + "\n" + seq + "\n")
    return primer_list


def split_list_average(list_in, number):
    result_list = []
    step = math.ceil(len(list_in) / number)
    for i in range(0, len(list_in), step):
        result_list.append(list_in[i:i + step])
    return result_list


if __name__ == '__main__':
    (options, args) = argsParse()
    print("INFO {} Start: Load file ...".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))
    product_len = options.len.split(",")
    print("off-targets length: {} - {}".format(product_len[0], product_len[1]))
    fasta = options.input
    e1 = time.time()
    term9_fasta_tmp = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".9.fa")
    # term9_fasta = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "9.fa"
    get_term9(fasta, term9_fasta_tmp)
    term9_fasta = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".expand.fa")
    primer_info = fa_dege_trans(term9_fasta_tmp, term9_fasta)
    path_2_ref = options.ref
    sam_file = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".sam")
    sam_for_file = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".for.sam")
    sam_rev_file = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".rev.sam")
    # sam_file = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "sam"
    # sam_for_file = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "for.sam"
    # sam_rev_file = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "rev.sam"
    if sam_for_file.exists() and sam_rev_file.exists():
        pass
    else:
        bowtie_map(term9_fasta, path_2_ref, sam_file, sam_for_file, sam_rev_file, options.proc)
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
    gene_list = split_list_average(target_gene, options.proc)
    # pool.apply_async(PCR_prediction, args=(a,), kwds={b: value})
    prediction = multiprocessing.Manager().list()
    pool2 = multiprocessing.Pool()
    for i in gene_list:
        pool2.apply_async(PCR_prediction, args=(i, forward_dict, reverse_dict, prediction))
    pool2.close()
    pool2.join()
    print("Number of predicted PCR products: {}.".format(len(prediction)))
    prediction = pd.DataFrame(list(prediction))
    with open(options.out, "w") as f:
        prediction.to_csv(f, index=False, sep="\t")
    ########################################################################################
    primer_number_out = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(
        ".predicted_product_number")
    primer_list = ["Primer_F", "Primer_F"]
    pieces = []
    for col in primer_list:
        tmp_series = prediction[col].value_counts()
        tmp_series.name = col
        pieces.append(tmp_series)
    df_value_counts = pd.concat(pieces, axis=1)
    df_value_counts.fillna(0, inplace=True)
    df_value_counts["rowSum"] = df_value_counts.apply(lambda x: x.sum(), axis=1)
    df_value_counts.sort_values(["rowSum"], ascending=False)
    with open(primer_number_out,"w") as f:
        df_value_counts.to_csv(f, index=True, sep="\t")
    print("INFO {} Done...".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
