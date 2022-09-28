import os
import argparse
import time
from math import log10
from itertools import product
from multiprocessing import Manager
from collections import defaultdict

from concurrent.futures import ProcessPoolExecutor

# TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")
TRANS = str.maketrans("ATGC", "TACG")

degenerate_base = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"], "K": ["G", "T"],
                   "S": ["G", "C"], "W": ["A", "T"], "H": ["A", "T", "C"], "B": ["G", "T", "C"],
                   "V": ["G", "A", "C"], "D": ["G", "A", "T"], "N": ["A", "T", "G", "C"]}

score_table = {"A": 1, "G": 1.1, "C": 1.2, "T": 1.4, "R": 2.1, "Y": 2.6, "M": 2.2,
               "K": 2.5, "S": 2.3, "W": 2.4, "H": 3.6, "B": 3.7, "V": 3.3, "D": 3.5, "N": 4.7}

freedom_of_H_37_table = [[-0.7, -0.81, -0.65, -0.65],
                         [-0.67, -0.72, -0.8, -0.65],
                         [-0.69, -0.87, -0.72, -0.81],
                         [-0.61, -0.69, -0.67, -0.7]]

penalty_of_H_37_table = [[0.4, 0.575, 0.33, 0.73],
                         [0.23, 0.32, 0.17, 0.33],
                         [0.41, 0.45, 0.32, 0.575],
                         [0.33, 0.41, 0.23, 0.4]]

H_bonds_number = [[2, 2.5, 2.5, 2],
                  [2.5, 3, 3, 2.5],
                  [2.5, 3, 3, 2.5],
                  [2, 2.5, 2.5, 2]]

adjust = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}

adjust_initiation = {"A": 2.8, "T": 2.8, "C": 1.82, "G": 1.82}

adjust_terminal_TA = 0.4

symmetry_correction = 0.4

base2bit = {"A": 0, "C": 1, "G": 2, "T": 3}


def reversecomplement(seq):
    return seq.translate(TRANS)[::-1]


def Penalty_points(length, GC, d1, d2):
    return log10((2 ** length * 2 ** GC) / ((d1 + 0.1) * (d2 + 0.1)))


def parseArg():
    parser = argparse.ArgumentParser(
        description="For primer dimer check")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="input fasta primer file", metavar="<file>")
    parser.add_argument("-n", "--num", type=int, default=5,
                        help='number of cpu process, 5 by default', metavar="<int>")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help='output file', metavar="<file>")
    return parser.parse_args()


class Dimer(object):

    def __init__(self, primer_file="", outfile="", nproc=10):
        self.nproc = nproc
        self.primers_file = primer_file
        self.outfile = os.path.abspath(outfile)
        self.primers = self.parse_primers()
        self.resQ = Manager().Queue()

    def parse_primers(self):
        primer_dict = defaultdict(str)
        with open(self.primers_file,"r") as f:
            for i in f:
                if i.startswith(">"):
                    name = i.strip()
                else:
                    primer_dict[i.strip()] = name
        return primer_dict

    @staticmethod
    def degenerate_seq(primer):
        seq = []
        cs = ""
        for s in primer:
            if s not in degenerate_base:
                cs += s
            else:
                seq.append([cs + i for i in degenerate_base[s]])
                cs = ""
        if cs:
            seq.append([cs])
        return ("".join(i) for i in product(*seq))

    def current_end(self, primer, adaptor="", num=5, length=14):
        primer_extend = adaptor + primer
        end_seq = []
        for i in range(num, (num + length)):
            s = primer_extend[-i:]
            if s:
                end_seq.extend(self.degenerate_seq(s))
        return end_seq

    def deltaG(self, sequence):
        Delta_G_list = []
        for seq in self.degenerate_seq(sequence):
            Delta_G = 0
            for n in range(len(seq) - 1):
                i, j = base2bit[seq[n + 1]], base2bit[seq[n]]
                Delta_G += freedom_of_H_37_table[i][j] * H_bonds_number[i][j] + penalty_of_H_37_table[i][j]
            term5 = sequence[-2:]
            if term5 == "TA":
                Delta_G += adjust_initiation[seq[0]] + adjust_terminal_TA + symmetry_correction
            else:
                Delta_G += adjust_initiation[seq[0]] + symmetry_correction
            Delta_G_list.append(Delta_G)
        return round(max(Delta_G_list), 2)

    def dimer_check(self, primer):
        current_end_set = self.current_end(primer)
        current_end_list = sorted(list(current_end_set), key=lambda i: len(i), reverse=True)
        dimer = False
        for ps in self.primers:
            for end in current_end_list:
                for p in self.degenerate_seq(ps):
                    idx = p.find(reversecomplement(end))
                    if idx >= 0:
                        end_length = len(end)
                        end_GC = end.count("G") + end.count("C")
                        end_d1 = 0
                        end_d2 = len(p) - len(end) - idx
                        Loss = Penalty_points(
                            end_length, end_GC, end_d1, end_d2)
                        delta_G = self.deltaG(end)
                        if Loss > 3 or delta_G < -5:
                            line = (self.primers[primer], primer,
                                    end, delta_G, end_length, end_d1,
                                    end_GC, self.primers[ps],
                                    ps, end_d2, Loss
                                    )
                            self.resQ.put(line)
                            # The put method also has two optional parameters: blocked and timeout. If blocked is
                            # true (the default value) and timeout is positive, the method will block the time
                            # specified by timeout until there is space left in the queue. If the timeout occurs,
                            # a Queue.Full exception will be thrown. If blocked is false, but the queue is full,
                            # a Queue.Full exception will be thrown immediately.
                            dimer = True
                            if dimer:
                                break
                if dimer:
                    dimer = False
                    break
        self.resQ.put(None)  # Found a None, you can rest now

    #  The queue in multiprocessing cannot be used for pool process pool, but there is a manager in multiprocessing.
    #  Inter process communication in the pool uses the queue in the manager. Manager().Queue().
    #  Queue. qsize(): returns the number of messages contained in the current queue;
    #  Queue. Empty(): returns True if the queue is empty, otherwise False;
    #  Queue. full(): returns True if the queue is full, otherwise False;
    #  Queue. get(): get a message in the queue, and then remove it from the queue,
    #                which can pass the parameter timeout.
    #  Queue.get_Nowait(): equivalent to Queue. get (False).
    #                If the value cannot be obtained, an exception will be triggered: Empty;
    #  Queue. put(): add a value to the data sequence to transfer the parameter timeout duration.
    #  Queue.put_Nowait(): equivalent to Queue. get (False). When the queue is full, an error is reported: Full.

    #  Realize interprocess communication through pipe, and the performance of pipe is higher than that of queue.
    #  Pipe can only be used for communication between two processes.

    def run(self):
        p = ProcessPoolExecutor(self.nproc)
        for primer in self.primers.keys():
            p.submit(self.dimer_check, primer)
            #  This will submit all tasks to one place without blocking, and then each
            #  thread in the thread pool will fetch tasks.
        n = 0
        primer_id_sum = defaultdict(int)
        dimer_primer_id_sum = defaultdict(int)
        with open(self.outfile, "w") as fo:
            headers = ["Primer_ID", "Primer seq", "Primer end", "Delta G", "Primer end length", "End (distance 1)",
                       "End (GC)", "Dimer-primer_ID", "Dimer-primer seq", "End (distance 2)", "Loss"]
            fo.write("\t".join(headers) + "\n")
            while n < len(self.primers):
                res = self.resQ.get()
                # The get method can read and delete an element from the queue. Similarly, the get method has two
                # optional parameters: blocked and timeout. If blocked is true (the default value) and timeout is
                # positive, no element is retrieved during the waiting time, and a Queue is thrown Empty exception.
                # If blocked is false, there are two cases. If a value of Queue is available, return the value
                # immediately. Otherwise, if the queue is empty, throw a Queue.Empty exception immediately.
                if res is None:
                    n += 1
                    continue
                primer_id_sum[res[0]] += 1
                dimer_primer_id_sum[res[7]] += 1
                fo.write("\t".join(map(str, res)) + "\n")
            #  get results before shutdown. Synchronous call mode: call, wait for the return value, decouple, but slow.
        p.shutdown()
        # After I run the main, I don't care whether the sub thread is alive or dead. With this parameter, after all
        # the sub threads are executed, the main function is executed.
        # get results after shutdown. Asynchronous call mode: only call, unequal return values, coupling may exist,
        # but the speed is fast.

        with open(self.outfile + ".dimer_num", "w") as fo:
            fo.write("SeqName\tPrimer_ID\tDimer-primer_ID\tRowSum\n")
            for k in primer_id_sum.keys():
                p_id = primer_id_sum[k]
                d_id = dimer_primer_id_sum[k]
                RowSum = p_id + d_id
                fo.write("\t".join(map(str, [k, p_id, d_id, RowSum])) + "\n")


def main():
    args = parseArg()
    dimer_app = Dimer(primer_file=args.input,
                      outfile=args.output, nproc=args.num)
    dimer_app.run()


if __name__ == "__main__":
    e1 = time.time()
    main()
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))
