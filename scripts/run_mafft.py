#!/bin/python
"""mafft --auto in > out"""
#caution the chrom column of the position file must contain chr

__date__ = "2023-2-24"
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
import os
import optparse
import time
import re
from optparse import OptionParser
from collections import defaultdict


def argsParse():
    parser = OptionParser('Usage: %prog -i [input.fa] -o [output.fa] ')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file')
    parser.add_option('-o', '--out',
                      dest='out',
                      help='Out file')
    (options, args) = parser.parse_args()
    import sys
    from sys import argv
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


def run_mafft(Input, Output):
    info = Input.rstrip(".tfa").split("_")
    #print(info[-1])
    if info[-1] == "1" or info[-1] == 1:
        os.system("muscle -in {} -out {}".format(Input, Output))
    else:
        os.system("mafft --auto {} > {}".format(Input, Output))


if __name__ == "__main__":
    (options, args) = argsParse()
    In = options.input
    Out = options.out
    e1 = time.time()
    run_mafft(In, Out)
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
