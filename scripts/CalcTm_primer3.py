#!/bin/python
# Bug fixed.
# Sometimes, positons have no bases except "-", and columns in frequency array become [0,0,0,0],
# which is not proper for primer design, especially for Tm calculation,
# we repalce "-" in the start and stop region with the flanking nucloetides.
# we also set H and S to 0 in this situation to continue Tm calculation.


__date__ = "2023-6-6"
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


import time
from optparse import OptionParser
import sys
import primer3
import os
from pathlib import Path


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -f [format] -o [output]', version="%prog 0.0.1")

    parser.add_option('-i', '--input',
                      dest='input',
                      help='Primer file. One of the followed three types:\n '
                           'final_maxprimers_set.xls \n primer.fa \n primer_F,primer_R.')

    parser.add_option('-f', '--format',
                      dest='format',
                      default="fa",
                      type="str",
                      help='Format of primer file: fa or seq; default: fa. \n '
                           'fa: fasta format. \n'
                           'seq: sequence format. e.g. ATCCCG.')

    parser.add_option('-o', '--out',
                      dest='out',
                      default="primer_Tm.xls",
                      type="str",
                      help='Output. default: primer_Tm.xls.')



    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
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


class Cal_Tm(object):
    def __init__(self, primer_file="", output_file="", file_format="fa"):
        self.primers_file = primer_file
        self.output_file = output_file
        self.file_format = file_format


    def run(self):
        if self.file_format == "fa":
            with open(self.primers_file, "r") as f:
                with open(self.output_file, "w") as o:
                    for row in f:
                        if row.startswith("#") or row == "\n":
                            pass
                        elif row.startswith(">"):
                            primer_info = row.strip()
                        else:
                            primer = row.strip()
                            Tm = primer3.calcTm(primer)
                            o.write(primer_info + "\t" + primer + "\t" + str(Tm) + "\n")
        else:
            print("{}: {}".format(self.primers_file, primer3.calcTm(self.primers_file)))
            with open(self.output_file, "w") as o:
                Tm = primer3.calcTm(self.primers_file)
                o.write(self.primers_file + "\t" + str(Tm) + "\n")
def main():
    (options, args) = argsParse()
    results = Cal_Tm(primer_file=options.input, output_file=options.out, file_format=options.format)
    results.run()

if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
