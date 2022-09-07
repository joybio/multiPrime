#!/bin/python
# !/share/data3/yangjunbo/miniconda3/bin/python
__date__ = "2022-8-31"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import os
import sys
from optparse import OptionParser


def argsParse():
    parser = OptionParser('Usage: %prog -i [input] -o [output] -n [number]'
                          'Options: -s []')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input dir: list of accession number.')
    parser.add_option('-s', '--step',
                      dest='step',
                      default="5",
                      type="int",
                      help='Step length. Step between Primer_1_F : Primer_2_F')
    parser.add_option('-o', '--out',
                      dest='out',
                      help='Output directory')
    parser.add_option('-n', '--number',
                      dest='num',
                      help='Primer number of each cluster.')

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


def cal_number(Input, number):
    with open(Input, "r") as f:
        with open(number, "w") as out:
            for i in f:
                i = i.strip().split("\t")
                n = 1
                primer_number = 1
                primer = i[0].split("/")[-1].strip(".candidate.primers.txt")
                primer_out = options.out + "/" + primer + ".candidate.primers.fa"
                with open(primer_out, "w") as f:
                    while n < len(i):
                        start_stop = i[n + 4].split(":")
                        f.write(">" + primer + "_" + start_stop[0] + "_F\n" + i[n] + "\n>" + primer + "_" + start_stop[
                            1] + "_R\n" + i[n + 1] + "\n")
                        n += step
                    primer_number += 1
                out.write(primer + "\t" + str(primer_number) + "\n")


if __name__ == "__main__":
    (options, args) = argsParse()
    output = options.out
    number = options.num
    if os.path.isdir(output):
        print("Output directory is OK!!!\n")
    else:
        os.system("mkdir {}".format(output))
    step = options.step
    cal_number(options.input, number)

