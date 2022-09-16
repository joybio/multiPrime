#!/bin/python
"""
Output off-target PCR results. 
"""
__date__ = "2022-8-15"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import sys
from collections import defaultdict
import os
import threading
import time
import pandas as pd
from optparse import OptionParser
from pathlib import Path
import math
from functools import reduce
from operator import mul  #

# Path(path).parent, Path(path).name, Path(path).suffix, Path(path).stem, Path(path).iterdir(), Path(path).joinpath()
# Path(path).is_absolute(), Path(path).is_dir(), Path(path).is_file(), Path(path).exists(), Path(path).with_suffix


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -f [forward.sam] -r [reverse.sam] -p [50,2000]-o [output]')

    parser.add_option('-i', '--input',
                      dest='input',
                      help='input file: primer.fa.')

    parser.add_option('-r', '--ref',
                      dest='ref',
                      help='reference file: bowtie index.')

    parser.add_option('-p', '--product',
                      dest='product',
                      default="150,2000",
                      help='Length of PCR product, default: [150,2000].')

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
    with open(fa, "r") as f:
        with open(out, "w") as o:
            for i in f:
                if i.startswith(">"):
                    o.write(i)
                else:
                    i = i.strip()[-9:]
                    o.write(i + "\n")


def bowtie2_map(fa, ref_index, out, for_out, rev_out):
    # os.system("bowtie2 -p 20 -f -N 0 -a -x {} -f -U {} -S {}".format(ref_index, fa, out))
    os.system("bowtie -f -v 0 -n 1 -a -p 20 --best --strata {} {} -S {}".format(ref_index, fa, out))
    os.system("samtools view -F 16 {} > {}".format(out, for_out))
    os.system("samtools view -f 16 {} > {}".format(out, rev_out))


def build_dict(Input, Input_dict):
    with open(Input, "r") as f:
        for i in f:
            i = i.strip().split("\t")
            primer = i[0]
            gene = i[2]
            primer_match_start = int(i[3]) - 1
            Input_dict[gene].append([primer, primer_match_start])


def PCR_prediction(F_dict, R_dict):
    global prediction
    for gene in F_dict.keys():
        for primer_F in F_dict[gene]:
            position_start = primer_F[1]
            for primer_R in R_dict[gene]:
                position_stop = primer_R[1]
                distance = int(position_stop) - int(position_start) + 1
                if int(product_len[0]) < distance < int(product_len[1]):
                    offtarget = pd.DataFrame({"Primer_F": primer_F[0],
                                              "Primer_R": primer_R[0],
                                              "Product length": distance,
                                              "Chrom (or Genes)": gene,
                                              "Start": int(primer_F[1]) - 1,
                                              "Stop": int(primer_R[1]) - 1
                                              }, index=[0])
                    prediction = pd.concat([prediction, offtarget], axis=0, ignore_index=True)


def fa_dege_trans(fa, output):
    out = open(output, "w")
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
                        out.write(ID + "\n" + expand_seq[j] + "\n")
                else:
                    ID = Id + "_0"
                    out.write(ID + "\n" + expand_seq[0] + "\n")
    out.close()



if __name__ == '__main__':
    (options, args) = argsParse()
    print("INFO {} Start: Load file ...".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))
    product_len = options.product.split(",")
    print("off-targets length: {} - {}".format(product_len[0], product_len[1]))
    fasta = options.input
    term9_fasta_tmp = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".9.fa")
    # term9_fasta = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "9.fa"
    get_term9(fasta, term9_fasta_tmp)
    term9_fasta = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".expand.fa")
    fa_dege_trans(term9_fasta_tmp, term9_fasta)
    path_2_ref = options.ref
    sam_file = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".sam")
    sam_for_file = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".for.sam")
    sam_rev_file = Path(options.out).parent.joinpath(Path(options.input).stem).with_suffix(".rev.sam")
    # sam_file = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "sam"
    # sam_for_file = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "for.sam"
    # sam_rev_file = str(Path(options.out).parent) + "/" + os.path.basename(options.input).rstrip("fa") + "rev.sam"
    bowtie2_map(term9_fasta, path_2_ref, sam_file, sam_for_file, sam_rev_file)
    forward_dict = defaultdict(list)
    reverse_dict = defaultdict(list)
    t = [threading.Thread(target=build_dict, args=(sam_for_file, forward_dict)),
         threading.Thread(target=build_dict, args=(sam_rev_file, reverse_dict))]
    for t1 in t:
        t1.start()
    for t1 in t:
        t1.join()
    # print(forward_dict)
    prediction = pd.DataFrame(columns=["Primer_F", "Primer_R", "Product length", "Chrom (or Genes)", "Start", "Stop"])
    PCR_prediction(forward_dict, reverse_dict)
    # print(prediction)
    with open(options.out, "w") as f:
        prediction.to_csv(f, index=False, sep="\t")
    print("INFO {} Done...".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))
