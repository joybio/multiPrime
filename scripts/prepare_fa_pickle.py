#!/share/data3/yangjunbo/miniconda3/bin/python
#caution the chrom column of the position file must contain chr
__date__ = "2022-6-24"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import pickle
import optparse
from optparse import OptionParser
import re

parser = OptionParser('Usage: %prog -i [input] -o [output] \n return a binary dict, which can loaded by pickle')
parser.add_option('-i','--input',
		dest='input',
		help='Input file: Fasta. Caution: this fasta file must be formated, ID \n Sequence \n .......\n \
			make sure that you have run seq_format.py before, and use the output file as input.')
parser.add_option('-t','--head',
                dest='head',
		default="F",
                help='head information of the fasta, only use the title as the value of dict. default: F.key = NCBI ID/accession number. value = title')

parser.add_option('-o','--out',
		dest='out',
		help='Out file: Dictionary: key = NCBI ID/accession number. value = title \n sequence')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

data = open(options.input,"r")
out = open(options.out,"wb")
pickle_dict = {}
headinfo = options.head
for i in data:
	if i.startswith(">"):
		line = i
		i = i.lstrip(">").strip()
		key_list = []
		if re.search(">",i):
		#>FJ608688.1 Human parainfluenza virus 4a isolate 20 phosphoprotein gene, partial cds >FJ608691.1 Human parainfluenza virus 4a isolate 23 phosphoprotein gene, partial cds
			ID_list = i.split(">")
			for j in ID_list:
				j = j.split(" ")
				key_list.append(j[0])
		else:
			i = i.split(" ")
			key_list.append(i[0])
	else:
		if headinfo == "F":
			value = line + i
			#value = line
			for k in key_list:
				pickle_dict[k] = value
		else:
			value = line.strip()
			for k in key_list:
				pickle_dict[k] = value

		
pickle.dump(pickle_dict,out)
data.close()
out.close()


