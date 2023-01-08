#!/bin/python

__date__ = "2022-12-23"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

from optparse import OptionParser
from collections import defaultdict
import time
import sys

def argsParse():
    parser = OptionParser('Usage: %prog -i input -e experiment -o output')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: prediction of primer-dimer.')

    parser.add_option('-e', '--exp',
                      dest='exp',
                      help='Input file: expreiment with prediction primer set.')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='Output file: correlation.')
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Input file 1 must be specified !!!")
        sys.exit(1)
    elif options.exp is None:
        parser.print_help()
        print("Input file 2 must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("No output file provided !!!")
        sys.exit(1)
    return parser.parse_args()

def correlation(file1, file2, out):
	file1_loss_function_dict, file1_deltaG_dict, file1_num, file2_dict =defaultdict(int), defaultdict(int), defaultdict(int), defaultdict(int)
	with open(file1,"r") as f1:
		for i in f1:
			if i.startswith("Primer"):
				pass
			else:
				i = i.strip().split("\t")
				dimer = [i[0].lstrip(">"),i[7].lstrip(">")]
				dimer_sort = sorted(dimer)
				key = ' | '.join(dimer_sort)
				file1_loss_function_dict[key] += float(i[10])
				file1_deltaG_dict[key] += float(i[3])
				file1_num[key] += 1
	
	with open(file2,"r") as f2:
		for j in f2:
			j = j.strip().split("\t")
			dimer = [j[0],j[1]]
			dimer_sort = sorted(dimer)
			key = ' | '.join(dimer_sort)
			file2_dict[key] += int(j[2])

	with open(out,"w") as o:
		for i in file1_loss_function_dict.keys():
			if i in file2_dict.keys():
				o.write(i+"\t"+str(round(file1_loss_function_dict[i]/file1_num[i],2)) +"\t"+ \
					str(round(file1_deltaG_dict[i]/file1_num[i],2)) +"\t"+str(round(file2_dict[i],2)) +"\n")
			else:
				o.write(i+"\t"+str(round(file1_loss_function_dict[i]/file1_num[i],2)) +"\t"+ \
					str(round(file1_deltaG_dict[i]/file1_num[i],2)) +"\t0" +"\n")


def main():
	options, args = argsParse()
	correlation(options.input, options.exp, options.out)
if __name__ == "__main__":
	e1 = time.time()
	main()
	e2 = time.time()
	print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))


