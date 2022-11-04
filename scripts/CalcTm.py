# This CalcTm.py is similar with primer3
# reference:
# A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics
# Predicting DNA duplex stability from the base sequence
# Effects of Sodium Ions on DNA Duplex Oligomers: Improved Predictions of Melting Temperatures
# Predicting Stability of DNA Duplexes in Solutions Containing Magnesium and Monovalent Cations
# Oligonucleotide Melting Temperatures under PCR Conditions: Nearest-Neighbor Corrections for Mg21 , Deoxynucleotide
#       Triphosphate, andDimethyl Sulfoxide Concentrations with Comparison to Alternative Empirical Formulas

import math
import sys
from optparse import OptionParser


def argsParse():
    parser = OptionParser('Usage: %prog .')
    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input sequence: for example: "ATCTTCTCTC".')

    parser.add_option('-p', '--primer_conc',
                      dest='primer_conc',
                      default="50",
                      type="int",
                      help='primer concentration. Default: 50 nM')

    parser.add_option('-m', '--mono_conc',
                      dest='mono_conc',
                      default="50",
                      type="int",
                      help="monovalent concentration. Default: 50 nM.")

    parser.add_option('-d', '--diva_conc',
                      dest='diva_conc',
                      default="1.5",
                      type="float",
                      help="divalent concentration. Default: 1.5 mM.")

    parser.add_option('-n', '--dntp_conc',
                      dest='dntp_conc',
                      default="0.25",
                      type="float",
                      help="dntp concentration. Default: 0.25 mM.")

    parser.add_option('-o', '--out',
                      dest='out',
                      default="Tm.out",
                      help='Output file.')
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Input file must be specified !!!")
        sys.exit(1)
#    elif options.out is None:
#        parser.print_help()
#        print("No output file provided !!!")
#        sys.exit(1)
    return parser.parse_args()


# In duplex DNA there are 10 such unique doublets.
# These are 5 - 3 AA=TT; AG=CT; AC=GT; GA=TC; GG=CC; TG=CA, CG, GC, AT, and TA
# 1 cal = 4.184 J

# /* Table 1 (old parameters):
#  * See table 2 in the paper [Breslauer KJ, Frank R, Blöcker H and
#  * Marky LA (1986) "Predicting DNA duplex stability from the base
#  * sequence" Proc Natl Acad Sci 83:4746-50
#  * http://dx.doi.org/10.1073/pnas.83.11.3746]
#  */
# 1 M NaCl, 25 °C, and pH 7

# Htable1 = [[9.1, 5.8, 5.6, 6],
#            [6.5, 11, 11.1, 5.6],
#            [7.8, 11.9, 11, 5.8],
#            [8.6, 7.8, 6.5, 9.1]]
#
# Stable1 = [[24, 12.9, 13.5, 16.9],
#            [17.3, 26.6, 26.7, 13.5],
#            [20.8, 27.8, 26.6, 12.9],
#            [23.9, 20.8, 17.3, 24]]
#
# Gtable1 = [[1.9, 1.9, 1.6, 0.9],
#            [1.3, 3.1, 3.1, 1.6],
#            [1.6, 3.6, 3.1, 1.9],
#            [1.5, 1.6, 1.3, 1.9]]

# /* Table 2, new parameters:
#  * Tables of nearest-neighbor thermodynamics for DNA bases, from the
#  * paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
#  * and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
#  * Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]
#  */

TRANS = str.maketrans("ATCG", "TAGC")


def complement(seq):
    return seq.translate(TRANS)[::-1]


# 37°C and 1 M NaCl
Htable2 = [[-7.9, -8.5, -8.2, -7.2],
           [-8.4, -8, -9.8, -8.2],
           [-7.8, -10.6, -8, -8.5],
           [-7.2, -7.8, -8.4, -7.9]]

Stable2 = [[-22.2, -22.7, -22.2, -21.3],
           [-22.4, -19.9, -24.4, -22.2],
           [-21, -27.2, -19.9, -22.7],
           [-20.4, -21, -22.4, -22.2]]

Gtable2 = [[-1, -1.45, -1.3, -0.58],
           [-1.44, -1.84, -2.24, -1.3],
           [-1.28, -2.17, -1.84, -1.45],
           [-0.88, -1.28, -1.44, -1]]

base2bit = {"A": 0, "C": 1, "G": 2, "T": 3}

H_adjust_initiation = {"A": 2.3, "T": 2.3, "C": 0.1, "G": 0.1}
S_adjust_initiation = {"A": 4.1, "T": 4.1, "C": -2.8, "G": -2.8}
G_adjust_initiation = {"A": 1.03, "T": 1.03, "C": 0.98, "G": 0.98}
H_symmetry_correction = 0
S_symmetry_correction = -1.4
G_symmetry_correction = 0.4
Kelvin = 273.15
crossover_point = 0.22


def symmetry(seq):
    if len(seq) % 2 == 1:
        return False
    else:
        F = seq[:int(len(seq) / 2)]
        R = complement(seq[int(len(seq) / 2):][::-1])
        if F == R:
            return True
        else:
            return False


def Calc_deltaH_deltaS(seq):
    Delta_H = 0
    Delta_S = 0
    for n in range(len(seq) - 1):
        i, j = base2bit[seq[n + 1]], base2bit[seq[n]]
        Delta_H += Htable2[i][j]
        Delta_S += Stable2[i][j]
    Delta_H += H_adjust_initiation[seq[0]] + H_adjust_initiation[seq[-1]]
    Delta_S += S_adjust_initiation[seq[0]] + S_adjust_initiation[seq[-1]]
    if symmetry(seq):
        Delta_S += S_symmetry_correction
    return Delta_H * 1000, Delta_S


# salt_adjust = math.log(Tm_Na_adjust / 1000.0, math.e)
#
#
# def S_adjust(seq):
#     n = len(seq) - 1
#     # S_Na_adjust = 0.847 * n * salt_adjust
#     # Oligonucleotide Melting Temperatures under PCR Conditions: Nearest-Neighbor Corrections for
#     # Mg2+ , Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide Concentrations with
#     # Comparison to Alternative Empirical Formulas
#     S_Na_adjust = 0.368 * n * salt_adjust
#     # A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics
#     return S_Na_adjust
# where n is the total number of phosphates in the duplex divided by 2,
# This is equal to the oligonucleotide length minus 1.

def GC_fraction(seq):
    return round((list(seq).count("G") + list(seq).count("C")) / len(list(seq)), 3)


# Tm is calculated from the predicted deltaH and deltaS
# Tm = deltaH/(deltaS + R*lnCT)
# R is the gas constant (1.987 cal/Kmol)
# CT: the total oligonucleotide strand concentration
# For non-self-complementary molecules, CT is replaced by CT/4
# if the strands are in equal concentration or by (CA - CB)/2
# if the strands are at different concentrations, where CA and CB are the concentrations of the more concentrated
# and less concentrated strands, respectively
# def Calc_Tm_v1(seq):
#     delta_H, delta_S = Calc_deltaH_deltaS(seq)
#     # print(delta_H / (delta_S + 1.987 * math.log(primer_concentration / pow(10, 9),math.e)))
#     Tm = delta_H / (delta_S + S_adjust(seq) + 1.987 * math.log(primer_concentration / 4 * pow(10, 9), math.e)) \
#          + 0.345 * GC_fraction(seq) + salt_adjust * (17.0 - 0.135 * GC_fraction(seq)) - 550 / len(seq) - Kelvin
#     return Tm


# different salt corrections for monovalent (Owczarzy et al.,2004) and divalent cations (Owczarzy et al.,2008)
def Calc_Tm_v2(seq):
    delta_H, delta_S = Calc_deltaH_deltaS(seq)
    # Note that the concentrations in the following Eq is mmol/L, In all other equations,concentration are mol/L
    # Monovalent cations are typically present as K+ and Tris+ in PCR buffer,
    # K+ is similar to Na+ in regard to duplex stabilization
    # if Di_concentration > dNTP_concentration:
    #     Tm_Na_adjust = Mo_concentration + 120 * math.sqrt(Di_concentration - dNTP_concentration)
    # else:
    #     Tm_Na_adjust = Mo_concentration
    Tm_Na_adjust = Mo_concentration

    if dNTP_concentration >= Di_concentration:
        free_divalent = 0.00000000001
    else:
        free_divalent = (Di_concentration - dNTP_concentration) / 1000.0
    R_div_monov_ratio = (math.sqrt(free_divalent)) / (Mo_concentration / 1000)

    if R_div_monov_ratio < crossover_point:
        # use only monovalent salt correction, [equation 22] (Owczarzy et al., 2004)
        correction = (((4.29 * GC_fraction(seq)) - 3.95) * pow(10, -5) * math.log(Tm_Na_adjust / 1000.0, math.e)) \
                     + (9.40 * pow(10, -6) * (pow(math.log(Tm_Na_adjust / 1000.0, math.e), 2)))
    else:
        # magnesium effects are dominant, [equation 16] (Owczarzy et al., 2008) is used
        # Table 2
        a = 3.92 * pow(10, -5)
        b = - 9.11 * pow(10, -6)
        c = 6.26 * pow(10, -5)
        d = 1.42 * pow(10, -5)
        e = - 4.82 * pow(10, -4)
        f = 5.25 * pow(10, -4)
        g = 8.31 * pow(10, -5)
        if R_div_monov_ratio < 6.0:
            a = 3.92 * pow(10, -5) * (
                    0.843 - (0.352 * math.sqrt(Tm_Na_adjust / 1000.0) * math.log(Tm_Na_adjust / 1000.0, math.e)))
            d = 1.42 * pow(10, -5) * (
                    1.279 - 4.03 * pow(10, -3) * math.log(Tm_Na_adjust / 1000.0, math.e) - 8.03 * pow(10, -3) * pow(
                math.log(Tm_Na_adjust / 1000.0, math.e), 2))
            g = 8.31 * pow(10, -5) * (
                    0.486 - 0.258 * math.log(Tm_Na_adjust / 1000.0, math.e) + 5.25 * pow(10, -3) * pow(
                math.log(Tm_Na_adjust / 1000.0, math.e), 3))
        # Eq 16
        correction = a + (b * math.log(free_divalent, math.e))
        + GC_fraction(seq) * (c + (d * math.log(free_divalent, math.e)))
        + (1 / (2 * (len(seq) - 1))) * (e + (f * math.log(free_divalent, math.e))
                                        + g * (pow((math.log(free_divalent, math.e)), 2)))

    if symmetry(seq):
        # Equation A
        Tm = round(1 / ((1 / (delta_H / (delta_S + 1.9872 * math.log(primer_concentration / (1 * pow(10, 9)), math.e))))
                        + correction) - Kelvin, 2)
    else:
        # Equation B
        Tm = round(1 / ((1 / (delta_H / (delta_S + 1.9872 * math.log(primer_concentration / (4 * pow(10, 9)), math.e))))
                        + correction) - Kelvin, 2)
    return Tm


if __name__ == "__main__":
    (options, args) = argsParse()
    Mo_concentration = options.mono_conc
    Di_concentration = options.diva_conc
    dNTP_concentration = options.dntp_conc
    primer_concentration = options.primer_conc
    melting_T = Calc_Tm_v2(options.input)
    print(melting_T)
    with open(options.out,"w") as f:
        f.write(options.input + "\t" + str(melting_T))
    f.close()
    sys.exit()
