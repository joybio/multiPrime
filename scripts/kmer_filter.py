#!/share/data3/yangjunbo/miniconda3/bin/python

from optparse import OptionParser
import sys
import re
import os
import os.path
from statistics import mean
import math
from functools import reduce
from operator import mul  #


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -o [output]')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: Input fasta.')
    parser.add_option('-g', '--gc',
                      dest='gc',
                      default="0.23,0.78",
                      help="Filter primers by GC content. Default [0.23,0.78].")
    parser.add_option('-o', '--out',
                      dest='out',
                      help='Output file: candidate kmer.')
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
    count = 3  # times for the packages install
    while count:
        try:
            import Bio  #
            print('Dependent package Biopython is OK.\nDependent module Bio is OK.')
            break
        except:
            print('Dependent package Biopython is not found!!! \n Start installing ....')
            os.system('pip install biopython')
            count -= 1
            continue
    return parser.parse_args()


########################## step1 screen #############################
#########################################################################################
# replace degenerate base to [A,C,G,T].


'''
def get_dict_key(d_dict, value):
	d_key = list(d_dict.keys())[list(d_dict.values()).index(value)]  # not for multi keys corresponding 1 value
	return d_key
'''
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


###########################################################

def GC_fraction(sequence):
    sequence_expand = dege_trans(sequence)
    GC_list = []
    for seq in sequence_expand:
        GC_list.append(round((list(seq).count("G") + list(seq).count("C")) / len(list(seq)), 3))
    GC_average = mean(GC_list)
    return GC_average


###########################################################
di_nucleotides = set()
bases = ["A", "C", "G", "T"]
for i in bases:
    single = i * 4
    di_nucleotides.add(single)
    for j in bases:
        if i != j:
            di = (i + j) * 4
            di_nucleotides.add(di)


def di_nucleotide(primer):
    Check = "False"
    primers = dege_trans(primer)
    for m in primers:
        for n in di_nucleotides:
            if re.search(n, m):
                Check = "True"
                break
            else:
                pass
        if Check == "True":
            break
    if Check == "True":
        return True
    else:
        return False


###########################################################
def GC_clamp(primer):
    end_5_nucleotide = primer[-5:]
    GC_percentage = GC_fraction(end_5_nucleotide)
    if GC_percentage > 0.6:
        return True
    else:
        return False


###########################################################

def pre_filter(Input, GC, Output):
    gc_content = GC
    with open(Input, "r") as In:
        with open(Output, "w") as Out:
            for seq in In:
                if seq.startswith(">"):
                    ID = seq
                else:
                    kmer = seq.strip()
                    GC_content = GC_fraction(kmer)
                    if di_nucleotide(kmer):
                        pass
                    elif GC_clamp(kmer):
                        pass
                    elif GC_content > float(gc_content[1]) or GC_content < float(gc_content[0]):
                        pass
                    else:
                        Out.write(ID + seq)


if __name__ == "__main__":
    (options, args) = argsParse()
    #### GC content ####
    GC = options.gc.split(",")
    print("GC content (candidate primer): {} - {}.\n".format(GC[0], GC[1]))
    pre_filter(options.input, GC, options.out)
