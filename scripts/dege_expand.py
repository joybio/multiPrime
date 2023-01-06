#!/bin/python
"""
expand degePrimers.
"""
__date__ = "2022-9-7"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import sys
from sys import argv
from optparse import OptionParser
from get_degePrimer import dege_trans
parser = OptionParser('Usage: %prog -i [input] -o [output]')
parser.add_option('-i', '--input',
                  dest='input',
                  help='Input file: degeprimer fa.')

parser.add_option('-o', '--out',
                  dest='out',
                  help='Output file: expand primers. e.g. expand degeprimer fa.')

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

def fa_dege_trans(fa, output):
    out = open(output, "w")
    with open(fa, "r") as dege_fa:
        for i in dege_fa:
            if i.startswith(">"):
                Id = i.strip()
            else:
                sequence = i.strip()
                expand_seq = dege_trans(sequence)
                for j,k in enumerate(expand_seq):
                        ID = Id + "_" + str(j)
                        out.write(ID + "\n" + k + "\n")
    out.close()

if __name__ == "__main__":
    tmp = options.input
    tmp_expand = options.out
    fa_dege_trans(tmp,tmp_expand)
        
       
