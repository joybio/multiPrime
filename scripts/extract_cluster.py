#!/bin/python

"""extract clusters from ch-hit results"""
# caution the chrom column of the position file must contain chr

__date__ = "2022-10-8"
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
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager
from optparse import OptionParser
import random
from collections import defaultdict
import re
import os

# TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")
from pathlib import Path

TRANS = str.maketrans("ATGC", "TACG")


def RC(seq):
    return seq.translate(TRANS)[::-1]


def argsParse():
    parser = OptionParser('Usage: %prog -i [total.fa] -c [cluster.clstr] \n'
                          'Options: -o [cluster.txt] -y [cluster.identities] -m [500] -d [Cluster.fa]',version = "%prog 0.0.2")
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: *complete.genome.fa')
    parser.add_option('-c', '--cluster',
                      dest='cluster',
                      help='cluster file: *.clstr')
    parser.add_option('-o', '--out',
                      dest='out',
                      default="cluster.txt",
                      type="str",
                      help='output file: clusters information. Default: cluster.txt.'
                           'It will be used as input file of snakemake pipeline')
    parser.add_option('-y', '--identity',
                      dest='identity',
                      default="cluster.identities.txt",
                      type="str",
                      help='output file: clusters information. Default: cluster.identities.txt.'
                           'It will be used as input file of snakemake pipeline')
    parser.add_option('-m', '--max',
                      dest='max',
                      default=1000,
                      type="int",
                      help='max sequence number in 1 cluster. Default: 1000.')
    parser.add_option('-d', '--dir',
                      dest='dir',
                      default="Cluster_fa",
                      type="str",
                      help='directory of output fasta: clusters information. Default: Cluster_fa.')

    parser.add_option('-p', '--proc',
                      dest='proc',
                      default="10",
                      type="int",
                      help='Number of process to launch.  default: 10.')

    (options, args) = parser.parse_args()
    import sys
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Input file must be specified !!!")
        sys.exit(1)
    elif options.cluster is None:
        parser.print_help()
        print("Cluster information (result of cd-hit) must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("No output file provided !!!")
        sys.exit(1)
    return parser.parse_args()


class Extract_Cluster(object):

    def __init__(self, Sequence_file="", Cluster_file="", Outfile="cluster.txt", Identity_file="cluster.identities.txt",
                 Seq_number=1000, Cluster_fa="Cluster_fa", nproc=10):
        self.seq = Sequence_file
        self.clstr = Cluster_file
        self.outfile = Outfile
        self.identity_file = Identity_file
        self.seq_number = Seq_number
        self.Cluster_fa = Cluster_fa
        self.nproc = nproc
        self.resQ = Manager().Queue()
        self.clstr_dict = self.parse_cluster()[0]
        self.clstr_rep = self.parse_cluster()[1]
        self.seq_dict = self.parse_seq()

    def parse_dir(self):
        Input_dir = Path(self.Cluster_fa)
        if Input_dir.exists():
            pass
        else:
            os.system("mkdir -p {}".format(Input_dir))

    def parse_cluster(self):
        cluster_rep = {}  # Representative sequence
        cluster_dict = defaultdict(list)
        cluster_identity = defaultdict(list)  # identity of each sequence (compare with representative sequence)
        with open(self.clstr, "r") as cluster:
            for c in cluster:
                if c.startswith(">"):
                    cluster_name = c.strip().replace(" ", "_")
                else:
                    c = c.strip().split(">")
                    value = c[1].split("...")
                    cluster_dict[cluster_name].append(">" + value[0])
                    identity = c[1].split(" ")
                    cluster_identity[cluster_name].append([value[0], identity[-1]])
                    if identity[-1] == "*":
                        cluster_rep[cluster_name] = ">" + value[0]
            cluster.close()
        with open(self.identity_file, "w") as identities:
            for i in cluster_identity.keys():
                for j in cluster_identity[i]:
                    if j[1] != "*":
                        identities.write(i.lstrip(">") + '\t' + j[0] + "\t" + j[1] + "\n")
            identities.close()
        identities.close()
        return [cluster_dict, cluster_rep]

    def parse_seq(self):
        seq_dict = {}
        with open(self.seq, "r") as data:
            for seq in data:
                if seq.startswith(">"):
                    if re.search(" ", seq):
                        seq = seq.split(" ")
                        seq_name = seq[0]
                    else:
                        seq_name = seq.strip()
                else:
                    seq_dict[seq_name] = seq
        data.close()
        return seq_dict

    def extract(self, ClusterID, number):
        current_cluster_file_id = ClusterID.lstrip(">") + "_" + str(len(self.clstr_dict[ClusterID])) + ".fa"
        current_cluster_file = Path(self.Cluster_fa).joinpath(current_cluster_file_id)
        with open(current_cluster_file, "w") as cf:
            for acc in self.clstr_dict[ClusterID]:
                cf.write(str(acc) + "\n" + self.seq_dict[acc])
        cf.close()
        top_current_cluster_file = Path(current_cluster_file).with_suffix(".tfa")
        current_cluster_seq = Path(current_cluster_file).with_suffix(".txt")
        with open(current_cluster_seq, "w") as cs:
            for i in self.clstr_dict[ClusterID]:
                cs.write(i.lstrip(">") + "\n")
        cs.close()
        if number == 0:
            os.system("cp {} {}".format(current_cluster_file, top_current_cluster_file))
        elif 0 < len(self.clstr_dict[ClusterID]) < int(number):
            with open(top_current_cluster_file, "w") as tmp:
                for acc_top in self.clstr_dict[ClusterID]:
                    tmp.write(str(acc_top) + "\n" + self.seq_dict[acc_top])
            tmp.close()
        else:
            self.clstr_dict[ClusterID].remove(self.clstr_rep[ClusterID])
            selected_seq = random.sample(self.clstr_dict[ClusterID], int(number) - 1)
            selected_seq.append(self.clstr_rep[ClusterID])
            with open(top_current_cluster_file, "w") as t:
                for seq_id in set(selected_seq):
                    t.write(str(seq_id) + "\n" + self.seq_dict[seq_id])
            t.close()

    def run(self):
        self.parse_dir()
        p = ProcessPoolExecutor(self.nproc)
        for ClusterID in self.clstr_dict.keys():
            with open(self.outfile, "w") as o:
                o.write("#Cluster_id\tNumber\n")
                for clusterID in self.clstr_dict.keys():
                    o.write(clusterID.lstrip(">") + "\t" + str(len(self.clstr_dict[clusterID])) + "\n")
            p.submit(self.extract, ClusterID, self.seq_number)
            #  This will submit all tasks to one place without blocking, and then each
            #  thread in the thread pool will fetch tasks.
        o.close()
        p.shutdown()
        # After I run the main, I don't care whether the sub thread is alive or dead. With this parameter, after all
        # the sub threads are executed, the main function is executed.
        # get results after shutdown. Asynchronous call mode: only call, unequal return values, coupling may exist,
        # but the speed is fast.


def main():
    (options, args) = argsParse()
    extract_cluster = Extract_Cluster(Sequence_file=options.input, Cluster_file=options.cluster, Outfile=options.out,
                                      Identity_file=options.identity, Seq_number=options.max, Cluster_fa=options.dir,
                                      nproc=options.proc)
    extract_cluster.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
