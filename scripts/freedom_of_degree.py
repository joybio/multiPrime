#
import pandas as pd

"""Calculate delta G of primer or kmer."""

__date__ = "2022-9-9"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

from optparse import OptionParser
import sys


def argsParse():
    parser = OptionParser('Usage: %prog -i input] -o [output] \n'
                          'Options: -g [gini] -p [position]')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file.')
    parser.add_option('-g', '--gini',
                      dest='gini',
                      default="unified",
                      type="str",
                      help='method used to calculate delta G. [unidfied] or [H_bonds]. Default:unified.\n'
                           'NN model: nearest neighbor model')
    parser.add_option('-p', '--position',
                      dest='position',
                      default="0",
                      type="int",
                      help='which column is sequence. Default: [0]')
    parser.add_option('-f', '--format',
                      dest='format',
                      help='Format of primer file: xls or fa or seq; default: xls. \n \
    				 fa: fasta format. \n seq: sequence format, e.g. ATGCTGATGCATCGT.')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='output file.')

    (options, args) = parser.parse_args()
    import sys
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


# Melting Temperature Tm(K)={ΔH/ ΔS + R ln(C)}, Or Melting Temperature Tm(oC) = {ΔH/ ΔS + R ln(C)} - 273.1
# ΔH (kcal/mole) : H is the Enthalpy. Enthalpy is the amount of heat energy possessed by substances. ΔH is the change
# in Enthalpy. In the above formula the ΔH is obtained by adding up all the di-nucleotide pairs enthalpy values of
# each nearest neighbor base pair. ΔS (kcal/mole) : S is the amount of disorder a system exhibits is called entropy.
# ΔS is change in Entropy. Here it is obtained by adding up all the di-nucleotide pairs entropy values of each
# nearest neighbor base pair. An additional salt correction is added as the Nearest Neighbor parameters were obtained
# from DNA melting studies conducted in 1M Na+ buffer and this is the default condition used for all calculations.
# ΔS (salt correction) = ΔS (1M NaCl )+ 0.368 x N x ln([Na+])
# Where N is the number of nucleotide pairs in the primer ( primer length -1).
# [Na+] is salt equivalent in mM.
# [Na+] calculation:
# [Na+] = Monovalent ion concentration +4 x free Mg2+.

# 1 cal = 4.184 J
# gas constant = 1.987 cal/Kmol) = 8.314J/(mol·K)
# Delta S = 24.85
# Delta H = T(i) * Delta S. T(i), where i are the NN dimers,Kelvin temperatures (273.15)
# Delta G = Delta H - T * (Delta S) || ΔG = ΔH – TΔS
# The dimer D G°37 values are then calculated with G = H-TS
Nearest_neighbor_stacking_energie_T_table = pd.DataFrame({"A": [63.88, 104.72, 76.5, 70.82],
                                                          "C": [68.7, 99.45, 88.74, 76.5],
                                                          "G": [89.45, 134.92, 99.45, 104.72],
                                                          "T": [56.31, 89.45, 68.7, 63.88]},
                                                         index=["A", "C", "G", "T"])
# AA: (63.88 + 273.15) * 24.85 - (60 + 273.15) * 24.85 =8375.2 - 8278.8 = 96.42 cal/mol = 0.096 kcal/mol
# Breslauer G°25
freedom_of_degree_25_table = pd.DataFrame({"A": [-1.9, -1.9, -1.6, -1.0],
                                           "C": [-1.3, -3.1, -3.1, -1.6],
                                           "G": [-1.6, -3.6, -3.1, -1.9],
                                           "T": [-1.5, -1.6, -1.3, -1.9]},
                                          index=["A", "C", "G", "T"])

# SantaLucia G°37  Unified. This work suggested that oligonucleotides with terminal 5' -T-A-3' base pairs should have
# a penalty of 1 0.4 kcal/mol but that no penalty should be given for terminal 5' -A-T-3' pair
freedom_of_degree_37_table_unified = pd.DataFrame({"A": [-1.00, -1.44, -1.28, -0.88],
                                                   "C": [-1.45, -1.84, -2.17, -1.28],
                                                   "G": [-1.30, -2.24, -1.84, -1.44],
                                                   "T": [-0.58, -1.30, -1.45, -1.00]},
                                                  index=["A", "C", "G", "T"])

# Using the current model with a unique mean pairing contribution per formed H-bond (−0.72 kcal·mol−1), we need
# to add the penalties for each step and add the number of formed H-bonds (note, half-contribution of the
# H-bonds from each base pair in each bp step)
freedom_of_H_37_table = pd.DataFrame({"A": [-0.70, -0.67, -0.69, -0.61],
                                      "C": [-0.81, -0.72, -0.87, -0.69],
                                      "G": [-0.65, -0.80, -0.72, -0.67],
                                      "T": [-0.65, -0.65, -0.81, -0.70]},
                                     index=["A", "C", "G", "T"])

penalty_of_H_37_table = pd.DataFrame({"A": [0.4, 0.23, 0.41, 0.33],
                                      "C": [0.575, 0.32, 0.45, 0.41],
                                      "G": [0.33, 0.17, 0.32, 0.23],
                                      "T": [0.73, 0.33, 0.575, 0.4]},
                                     index=["A", "C", "G", "T"])

H_bonds_number = pd.DataFrame({"A": [2, 2.5, 2.5, 2],
                               "C": [2.5, 3, 3, 2.5],
                               "G": [2.5, 3, 3, 2.5],
                               "T": [2, 2.5, 2.5, 2]},
                              index=["A", "C", "G", "T"])

# print(freedom_of_degree_25_table)
# print(penalty_of_H_37_table.loc["G"])
# print(H_bonds_number.loc["G"])


adjust = {"A": 0.98, "T": 0.98, "C": 1.03, "G": 1.03}


def delta_G(sequence, method):
    i = 0
    Delta_G = 0
    # G2_average = 0
    if method == "unified":
        while i < len(sequence) - 1:
            Delta_G += freedom_of_degree_37_table_unified.loc[sequence[i + 1], sequence[i]]
            i += 1
        start = sequence[0]
        stop = sequence[-1]
        Delta_G += adjust[start] + adjust[stop]
    elif method == "H_bonds":
        while i < len(sequence) - 1:
            Delta_G += (freedom_of_H_37_table.loc[sequence[i + 1], sequence[i]] * H_bonds_number.loc[
                sequence[i + 1], sequence[i]]) + \
                       penalty_of_H_37_table.loc[sequence[i + 1], sequence[i]]
            i += 1
        start = sequence[0]
        stop = sequence[-1]
        Delta_G += adjust[start] + adjust[stop]
    return round(Delta_G, 2)


if __name__ == "__main__":
    (options, args) = argsParse()
    if options.gini != "unified" and options.gini != "H_bonds":
        print("Method is wrong")
        sys.exit(1)
    else:
        if options.format == 'seq':
            sequence = options.input
            Delta_G = delta_G(sequence, options.gini)
            with open(options.out, "w") as output:
                output.write(sequence + "\t" + str(Delta_G) + "\n")
            print(Delta_G)
        else:
            f = open(options.input, "r")
            if options.format == 'xls':
                with open(options.out, "w") as output:
                    for i in f:
                        line = i.strip()
                        i = i.strip().split("\t")
                        sequence = i[options.position]
                        Delta_G = delta_G(sequence, options.gini)
                        output.write(line + "\t" + str(Delta_G) + "\n")
            elif options.format == 'fa':
                with open(options.out, "w") as output:
                    for i in f:
                        if i.startswith(">"):
                            output.write(i.strip() + "\t")
                        else:
                            sequence = i.strip()
                            Delta_G = delta_G(sequence, options.gini)
                            output.write(i + "\t" + str(Delta_G) + "\n")
            else:
                print("Check your input!!!")
            f.close()
