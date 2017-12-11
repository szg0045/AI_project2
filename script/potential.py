import argparse
import os
import subprocess
import string
from sklearn import *
import numpy
from scipy import spatial
from scipy.stats.stats import pearsonr
import math

# $amino_acids = "ACDEFGHIKLMNPQRSTVWY"; #order used in encoding for SVM or NN
aa_levitt = "GAVLIPDENQKRSTMCYWHF"  # aa order used in Levitt contact potential map
aa_braun = "GAVLIFYWMCPSTNQHKRDE"  # aa order used in Braun contact potential map


###########Contact Potential Database#######################################################################
# We use the following contact potentials:
# 1. Hinds and Levitt (JMB) (two printed copies)  [used]
# 2. Miyazawa and Jernigan(JMB) (two printed copies) [used]
# 3. Zhang and Kim (PNAS) (alpha/beta/coil contact energy) [not used, may added later]
# 4. Zhu and Braun (Protein Sci, beta-strand only) [not used, may added later]
###########################################################################################################

# Hinds and Levitt
#	Gly  Ala  Val  Leu  Ile  Pro  Asp  Glu  Asn  Gln  Lys  Arg  Ser  Thr  Met  Cys  Tyr  Trp  His  Phe
# Gly	.1  .7   .1   .1   0    .5   .4   .6   .1   0    .4   -0.1 .4   .2   -0.1 -0.1 -0.4 -0.7 0    -0.3
# Ala        .5   -0.3 -0.4 -0.4 .6   .3   .6   .3   0    1.0  .2   .5   0    -0.5 .3   -0.7 -0.8 0    -0.8
# Val             -1.1 -1.2 -1.2 0    .4   0    0    -0.4 0.1  -0.5 0    -0.3 -1.0 -0.5 -1.2 -1.6 -0.5 -1.5
# Leu                  -1.4 -1.4 -0.1 0    -0.1 -0.1 -0.6 0.1  -0.6 0    -0.3 -1.3 -0.8 -1.4 -1.7 -0.7 -1.6
# Ile                       -1.5 -0.1 0    -0.2 -0.1 -0.4 0    -0.7 -0.1 -0.6 -1.4 -0.8 -1.4 -1.8 -0.8 -1.7
# Pro                            .1   .1   .1   -0.1 -0.3 .6   -0.2 .2   0    -0.5 0    -1.0 -1.3 -0.4 -0.7
# Asp                                 0    0    -0.6 -0.3 -1.0 -1.4 -0.3 -0.3 0.1  0    -1.0 -0.6 -1.1 -0.3
# Glu                                      0.1  -0.6 -0.4 -1.1 -1.5 -0.2 -0.3 -0.3 0.1  -1.0 -0.8 -1.0 -0.5
# Asn                                           -0.7 -0.7 -0.3 -0.8 -0.1 -0.4 -0.3 0    -0.8 -0.8 -0.8 -0.6
# Gln                                                -0.5 -0.4 -0.9 0    -0.5 -0.6 -0.2 -1.1 -1.0 -0.5 -0.8
# Lys                                                     .7   .1   .1   0    -0.1 .5   -1.0 -0.8 0    -0.4
# Arg                                                          -0.9 -0.4 -0.6 -0.5 0    -1.4 -1.3 -1.0 -0.9
# Ser                                                               0    -0.2 -0.1 -0.1 -0.6 -0.6 -0.6 -0.4
# Thr                                                                    -0.5 -0.6 -0.3 -0.8 -0.9 -0.7 -0.7
# Met                                                                         -1.5 -0.8 -1.5 -2.0 -0.9 -1.9
# Cys                                                                              -2.7 -0.8 -1.3 -0.6 -1.2
# Tyr                                                                                   -1.6 -1.8 -1.5 -1.7
# Trp                                                                                        -2.2 -1.5 -2.0
# His                                                                                             -1.6 -1.2
# Phe                                                                                                  -2.0
#################################################################################################################

# get Levitt contact potential
def getLevittCP(aa1, aa2):
    # my ($aa1, $aa2) = @_;
    # my $id1 = index($aa_levitt, $aa1);
    id1 = aa_levitt.find(aa1)
    # my $id2 = index($aa_levitt, $aa2);
    id2 = aa_levitt.find(aa2)
    if (id1 < 0 or id1 > 19):
        return 0
    if (id2 < 0 or id2 > 19):
        return 0
    if (id1 > id2):
        tmp = id1
        id1 = id2
        id2 = tmp

    # my @cp_map = ();
    cp_map = []
    for i in range(0, 20):
        cp_map.append("")
    cp_map[
        0] = "0.1  0.7   0.1   0.1   0    0.5   0.4   0.6   0.1   0    0.4   -0.1 0.4   0.2   -0.1 -0.1 -0.4 -0.7 0    -0.3"
    cp_map[1] = "0.5   -0.3 -0.4 -0.4 0.6   0.3   0.6   0.3   0    1.0  0.2   0.5   0    -0.5 .3   -0.7 -0.8 0    -0.8";
    cp_map[2] = "-1.1 -1.2 -1.2 0    0.4   0    0    -0.4 0.1  -0.5 0    -0.3 -1.0 -0.5 -1.2 -1.6 -0.5 -1.5";
    cp_map[3] = "-1.4 -1.4 -0.1 0    -0.1 -0.1 -0.6 0.1  -0.6 0    -0.3 -1.3 -0.8 -1.4 -1.7 -0.7 -1.6";
    cp_map[4] = "-1.5 -0.1 0    -0.2 -0.1 -0.4 0    -0.7 -0.1 -0.6 -1.4 -0.8 -1.4 -1.8 -0.8 -1.7";
    cp_map[5] = "0.1   0.1   0.1   -0.1 -0.3 0.6   -0.2 0.2   0    -0.5 0    -1.0 -1.3 -0.4 -0.7";
    cp_map[6] = "0    0    -0.6 -0.3 -1.0 -1.4 -0.3 -0.3 0.1  0    -1.0 -0.6 -1.1 -0.3";
    cp_map[7] = "0.1  -0.6 -0.4 -1.1 -1.5 -0.2 -0.3 -0.3 0.1  -1.0 -0.8 -1.0 -0.5";
    cp_map[8] = "-0.7 -0.7 -0.3 -0.8 -0.1 -0.4 -0.3 0    -0.8 -0.8 -0.8 -0.6";
    cp_map[9] = "-0.5 -0.4 -0.9 0    -0.5 -0.6 -0.2 -1.1 -1.0 -0.5 -0.8";
    cp_map[10] = "0.7   0.1   0.1   0    -0.1 0.5   -1.0 -0.8 0    -0.4";
    cp_map[11] = "-0.9 -0.4 -0.6 -0.5 0    -1.4 -1.3 -1.0 -0.9";
    cp_map[12] = "0    -0.2 -0.1 -0.1 -0.6 -0.6 -0.6 -0.4";
    cp_map[13] = "-0.5 -0.6 -0.3 -0.8 -0.9 -0.7 -0.7";
    cp_map[14] = "-1.5 -0.8 -1.5 -2.0 -0.9 -1.9";
    cp_map[15] = "-2.7 -0.8 -1.3 -0.6 -1.2";
    cp_map[16] = "-1.6 -1.8 -1.5 -1.7";
    cp_map[17] = "-2.2 -1.5 -2.0";
    cp_map[18] = "-1.6 -1.2";
    cp_map[19] = "-2.0";
    # my @energy = split(/\s+/, $cp_map[$id1]);
    energy = cp_map[id1].split()
    # my $cp = $energy[$id2 - $id1];
    cp = energy[id2 - id1]
    return float(cp)


# Miyazawa & Jernigan
#	Gly  Ala  Val  Leu  Ile  Pro  Asp  Glu  Asn  Gln  Lys  Arg  Ser  Thr  Met  Cys  Tyr  Trp  His  Phe
# Gly    -2.1 -2.2 -3.0 -2.5 -2.7 -1.8 -1.9 -1.3 -2.4 -2.0 -1.9 -2.2 -1.9 -2.4 -2.8 -3.0 -2.8 -3.1 -2.1 -2.6
# Ala         -2.9 -4.1 -3.7 -3.9 -2.3 -2.6 -1.9 -2.8 -2.7 -1.9 -2.4 -2.4 -3.2 -3.8 -3.1 -3.7 -3.8 -2.6 -3.7
# Val              -5.1 -4.7 -5.0 -3.2 -2.8 -2.8 -3.5 -3.4 -3.1 -3.5 -3.1 -3.9 -4.7 -4.3 -4.5 -4.8 -3.5 -4.6
# Leu                   -4.3 -4.6 -2.8 -2.7 -2.4 -3.0 -3.1 -2.5 -3.0 -2.6 -3.3 -4.4 -4.0 -4.1 -4.4 -3.1 -4.2
# Ile                        -4.9 -3.0 -2.9 -2.7 -3.2 -3.2 -2.9 -3.3 -3.0 -3.8 -4.7 -4.2 -4.4 -4.7 -3.4 -4.5
# Pro                             -2.2 -2.3 -1.8 -2.7 -2.5 -1.8 -2.4 -2.2 -2.8 -3.4 -2.8 -3.5 -3.7 -2.6 -3.1
# Asp                                  -2.6 -2.0 -3.3 -2.7 -3.5 -3.7 -2.8 -3.2 -2.7 -2.9 -3.5 -3.2 -3.3 -2.7
# Glu                                       -1.5 -2.8 -2.2 -3.2 -3.2 -2.3 -2.7 -2.7 -2.4 -3.0 -2.8 -2.8 -2.4
# Asn                                            -3.6 -3.2 -3.0 -3.2 -2.8 -3.4 -3.3 -3.2 -3.5 -3.5 -3.2 -3.2
# Gln                                                 -2.6 -2.7 -2.9 -2.3 -3.1 -3.3 -2.9 -3.4 -3.3 -2.5 -3.0
# Lys                                                      -1.7 -2.0 -2.3 -2.8 -2.9 -2.4 -3.5 -3.3 -2.2 -2.8
# Arg                                                           -2.9 -2.6 -3.1 -3.1 -2.6 -3.6 -3.5 -2.9 -3.0
# Ser                                                                -2.4 -2.9 -3.0 -3.0 -3.1 -3.1 -2.7 -2.8
# Thr                                                                     -3.5 -3.8 -3.5 -3.6 -3.7 -3.2 -3.4
# Met                                                                          -4.8 -4.2 -4.4 -4.8 -3.5 -4.6
# Cys                                                                               -6.1 -3.8 -4.3 -3.3 -4.1
# Tyr                                                                                    -4.1 -4.3 -3.7 -4.1
# Trp                                                                                         -4.7 -3.6 -4.4
# His                                                                                              -3.5 -3.3
# Phe                                                                                                   -4.3
##################################################################################################################
# get Jernigan contact potential
def getJerniganCP(aa1, aa2):
    # my ($aa1, $aa2) = @_;
    # my $id1 = index($aa_levitt, $aa1);
    id1 = aa_levitt.find(aa1)
    # my $id2 = index($aa_levitt, $aa2);
    id2 = aa_levitt.find(aa2)
    if (id1 < 0 or id1 > 19):
        return 0
    if (id2 < 0 or id2 > 19):
        return 0
    if (id1 > id2):
        tmp = id1
        id1 = id2
        id2 = tmp

    # my @cp_map = ();
    cp_map = []
    for i in range(0, 20):
        cp_map.append("")
    cp_map[0] = "-2.1 -2.2 -3.0 -2.5 -2.7 -1.8 -1.9 -1.3 -2.4 -2.0 -1.9 -2.2 -1.9 -2.4 -2.8 -3.0 -2.8 -3.1 -2.1 -2.6";
    cp_map[1] = "-2.9 -4.1 -3.7 -3.9 -2.3 -2.6 -1.9 -2.8 -2.7 -1.9 -2.4 -2.4 -3.2 -3.8 -3.1 -3.7 -3.8 -2.6 -3.7";
    cp_map[2] = "-5.1 -4.7 -5.0 -3.2 -2.8 -2.8 -3.5 -3.4 -3.1 -3.5 -3.1 -3.9 -4.7 -4.3 -4.5 -4.8 -3.5 -4.6";
    cp_map[3] = "-4.3 -4.6 -2.8 -2.7 -2.4 -3.0 -3.1 -2.5 -3.0 -2.6 -3.3 -4.4 -4.0 -4.1 -4.4 -3.1 -4.2";
    cp_map[4] = "-4.9 -3.0 -2.9 -2.7 -3.2 -3.2 -2.9 -3.3 -3.0 -3.8 -4.7 -4.2 -4.4 -4.7 -3.4 -4.5";
    cp_map[5] = "-2.2 -2.3 -1.8 -2.7 -2.5 -1.8 -2.4 -2.2 -2.8 -3.4 -2.8 -3.5 -3.7 -2.6 -3.1";
    cp_map[6] = "-2.6 -2.0 -3.3 -2.7 -3.5 -3.7 -2.8 -3.2 -2.7 -2.9 -3.5 -3.2 -3.3 -2.7";
    cp_map[7] = "-1.5 -2.8 -2.2 -3.2 -3.2 -2.3 -2.7 -2.7 -2.4 -3.0 -2.8 -2.8 -2.4";
    cp_map[8] = "-3.6 -3.2 -3.0 -3.2 -2.8 -3.4 -3.3 -3.2 -3.5 -3.5 -3.2 -3.2";
    cp_map[9] = "-2.6 -2.7 -2.9 -2.3 -3.1 -3.3 -2.9 -3.4 -3.3 -2.5 -3.0";
    cp_map[10] = "-1.7 -2.0 -2.3 -2.8 -2.9 -2.4 -3.5 -3.3 -2.2 -2.8";
    cp_map[11] = "-2.9 -2.6 -3.1 -3.1 -2.6 -3.6 -3.5 -2.9 -3.0";
    cp_map[12] = "-2.4 -2.9 -3.0 -3.0 -3.1 -3.1 -2.7 -2.8";
    cp_map[13] = "-3.5 -3.8 -3.5 -3.6 -3.7 -3.2 -3.4";
    cp_map[14] = "-4.8 -4.2 -4.4 -4.8 -3.5 -4.6";
    cp_map[15] = "-6.1 -3.8 -4.3 -3.3 -4.1";
    cp_map[16] = "-4.1 -4.3 -3.7 -4.1";
    cp_map[17] = "-4.7 -3.6 -4.4";
    cp_map[18] = "-3.5 -3.3";
    cp_map[19] = "-4.3";

    # my @energy = split(/\s+/, $cp_map[$id1]);
    energy = cp_map[id1].split()
    cp = energy[id2 - id1]
    return float(cp)


# Zhu-Braun
# G -0.29
# A -0.14 -0.18
# V -0.10 -0.15 -0.48
# L -0.04 -0.24 -0.29 -0.43
# I 0.27  -0.25 -0.31 -0.45 -0.48
# F -0.09 -0.16 -0.31 -0.28 -0.05 -0.50
# Y -0.21 -0.18 0.00  -0.10 -0.34 -0.27 -0.11
# W -0.34 -0.01 0.18  -0.18 -0.28 0.16  -0.30 -0.53
# M 0.25  -0.02 -0.02 -0.32 0.21  -0.36 0.01  -0.73 -0.75
# C -0.42 0.08  0.08  0.36  -0.16 -0.28 0.69  -0.74 0.27  -1.77
# P 0.06  0.28  0.76  0.30  0.99  0.65  -0.02 0.70  -0.78 0.31  -0.78
# S 0.04  0.38  0.18  0.30  0.57  0.15  -0.03 0.44  0.00  0.12  0.21  -0.68
# T 0.28  0.06  0.19  0.57  0.34  0.25  0.23  0.74  0.43  0.28  0.04  -0.23 -0.58
# N 0.49  -0.04 0.48  0.25  1.45  0.12  -0.14 0.46  -0.52 0.07  0.59  -0.21 -0.06 -0.45
# Q 0.54  0.35  0.41  0.35  0.44  -0.04 -0.06 -0.09 0.07  0.39  0.73  0.19  -0.31 0.20  -0.17
# H -0.09 0.44  0.37  0.10  0.24  0.25  0.33  -0.34 1.07  -0.45 -0.21 -0.13 -0.22 -0.56 0.28  -0.15
# K 0.56  0.28  0.53  0.37  -0.00 0.75  -0.00 0.02  0.44  0.68  0.26  -0.05 -0.26 -0.27 0.05  0.57  0.21
# R 0.40  0.59  0.43  0.37  0.05  0.31  0.03  -0.20 0.53  0.92  0.34  0.24  -0.31 -0.00 0.56  -0.11 0.58  -0.03
# D -0.26 0.24  0.51  0.80  0.26  0.33  0.61  0.74  0.21  0.53  0.87  -0.03 0.32  -0.43 -0.03 -0.61 -0.43 -0.79 0.11
# E 0.21  0.53  0.37  0.51  0.53  0.38  0.25  1.37  0.44  0.17  0.41  0.10  -0.27 0.76  -0.20 -0.14 -1.12 -0.85 0.86  0.58

def getBraunCP(aa1, aa2):
    # my ($aa1, $aa2) = @_;
    # my $id1 = index($aa_braun, $aa1);
    id1 = aa_braun.find(aa1)
    # my $id2 = index($aa_braun, $aa2);
    id2 = aa_braun.find(aa2)
    if (id1 < 0 or id1 > 19):
        return 0
    if (id2 < 0 or id2 > 19):
        return 0
    if (id1 < id2):
        tmp = id1
        id1 = id2
        id2 = tmp

    # my @cp_map = ();
    cp_map = []
    for i in range(0, 20):
        cp_map.append("")
    cp_map[0] = "-0.29";
    cp_map[1] = "-0.14 -0.18";
    cp_map[2] = "-0.10 -0.15 -0.48";
    cp_map[3] = "-0.04 -0.24 -0.29 -0.43";
    cp_map[4] = "0.27  -0.25 -0.31 -0.45 -0.48";
    cp_map[5] = "-0.09 -0.16 -0.31 -0.28 -0.05 -0.50";
    cp_map[6] = "-0.21 -0.18 0.00  -0.10 -0.34 -0.27 -0.11";
    cp_map[7] = "-0.34 -0.01 0.18  -0.18 -0.28 0.16  -0.30 -0.53";
    cp_map[8] = "0.25  -0.02 -0.02 -0.32 0.21  -0.36 0.01  -0.73 -0.75";
    cp_map[9] = "-0.42 0.08  0.08  0.36  -0.16 -0.28 0.69  -0.74 0.27  -1.77";
    cp_map[10] = "0.06  0.28  0.76  0.30  0.99  0.65  -0.02 0.70  -0.78 0.31  -0.78";
    cp_map[11] = "0.04  0.38  0.18  0.30  0.57  0.15  -0.03 0.44  0.00  0.12  0.21  -0.68";
    cp_map[12] = "0.28  0.06  0.19  0.57  0.34  0.25  0.23  0.74  0.43  0.28  0.04  -0.23 -0.58";
    cp_map[13] = "0.49  -0.04 0.48  0.25  1.45  0.12  -0.14 0.46  -0.52 0.07  0.59  -0.21 -0.06 -0.45";
    cp_map[14] = "0.54  0.35  0.41  0.35  0.44  -0.04 -0.06 -0.09 0.07  0.39  0.73  0.19  -0.31 0.20  -0.17";
    cp_map[15] = "-0.09 0.44  0.37  0.10  0.24  0.25  0.33  -0.34 1.07  -0.45 -0.21 -0.13 -0.22 -0.56 0.28  -0.15";
    cp_map[16] = "0.56  0.28  0.53  0.37  -0.00 0.75  -0.00 0.02  0.44  0.68  0.26  -0.05 -0.26 -0.27 0.05  0.57  0.21";
    cp_map[
        17] = "0.40  0.59  0.43  0.37  0.05  0.31  0.03  -0.20 0.53  0.92  0.34  0.24  -0.31 -0.00 0.56  -0.11 0.58  -0.03";
    cp_map[
        18] = "-0.26 0.24  0.51  0.80  0.26  0.33  0.61  0.74  0.21  0.53  0.87  -0.03 0.32  -0.43 -0.03 -0.61 -0.43 -0.79 0.11";
    cp_map[
        19] = "0.21  0.53  0.37  0.51  0.53  0.38  0.25  1.37  0.44  0.17  0.41  0.10  -0.27 0.76  -0.20 -0.14 -1.12 -0.85 0.86  0.58";

    # my @energy = split(/\s+/, $cp_map[$id1]);
    energy = cp_map[id1].split()
    # my $cp = $energy[$id2];
    cp = energy[id2]
    return float(cp)


# }

# test code###############
test = 0
if (test == 1):
    print "test levitt\n"
    print getLevittCP("A", "D"), "\n"
    print getLevittCP("P", "W"), "\n"
    print getLevittCP("R", "S"), "\n"
    print getLevittCP("K", "K"), "\n"
    print getLevittCP("L", "Y")

    print "test Jernigan \n"
    print getJerniganCP("A", "D"), "\n"
    print getJerniganCP("P", "W"), "\n"
    print getJerniganCP("R", "S"), "\n"
    print getJerniganCP("K", "K"), "\n"
    print getJerniganCP("L", "Y")

    print "test Braun \n"
    print getBraunCP("A", "D"), "\n"
    print getBraunCP("P", "W"), "\n"
    print getBraunCP("R", "S"), "\n"
    print getBraunCP("K", "K"), "\n"
    print getBraunCP("L", "Y")


def round(input_value):
    if isinstance(input_value, list):
        # print "list found: ", input_value
        value = input_value[0]
    else:
        # print type(input_value)
        # print input_value
        value = input_value
    value *= 100
    value = int(value + 0.5)
    value /= 100
    return value


# call: $value = &cosine(\@x,\@y)
# size of x, y array is 21
def cosine(x, y):
    # my ($x, $y) = @_;
    res = 0
    x_ave = 0
    y_ave = 0
    xy = 0
    size = len(x)
    if (size != len(y)):
        raise ValueError("size of two profiles doesn't equal in cosine.\n")
    x_len = 0
    y_len = 0
    for i in range(0, size):
        x_ave += x[i]
        y_ave += y[i]
        xy += (x[i] * y[i])
    x_ave /= size
    y_ave /= size
    for i in range(0, size):
        # $x_len += ($x->[$i] - $x_ave) * ($x->[$i] - $x_ave);
        # $y_len += ($y->[$i] - $y_ave) * ($y->[$i] - $y_ave);
        x_len += (x[i] * x[i])
        y_len += (y[i] * y[i])
    x_len = math.sqrt(x_len)
    y_len = math.sqrt(y_len)
    res = xy / (x_len * y_len)
    if (res < 0):
        res = 0
    if (res > 1):
        res = 1
    return res


def correlation(x, y):
    # my ($x, $y) = @_;
    res = 0
    x_ave = 0
    y_ave = 0
    xy = 0
    size = len(x)
    if (size != len(y)):
        raise ValueError("size of two profiles doesn't equal in correlation.\n")
    x_len = 0
    y_len = 0
    for i in range(0, size):
        x_ave += x[i]
        y_ave += y[i]
    x_ave /= size
    y_ave /= size
    for i in range(0, size):
        x_len += (x[i] - x_ave) * (x[i] - x_ave)
        y_len += (y[i] - y_ave) * (y[i] - y_ave)
        xy += (x[i] - x_ave) * (y[i] - y_ave)
    x_len = math.sqrt(x_len)
    y_len = math.sqrt(y_len)

    res = xy / (x_len * y_len)
    if (res < -1):
        res = -1
    if (res > 1):
        res = 1
    return res


# entropy
def entropy(x):
    # my $x = $_[0];
    ent = 0
    for i in range(0, len(x)):
        prob = x[i]
        if (prob > 0):
            ent -= (prob * math.log(prob))
    return ent


# mututal information
# three parameters: profile1(21), profile2(21), and joint dist(21*21)
# consider gap as one extra
def mutual(x, y, xy):
    import math
    # I(X:Y) = Sum_xy ( P(x,y)log(P(x,y)/p(x)p(y) )
    #:P(x,y): estimated from each pair frequence (21*21 combintions)
    # P(x), P(y): taken from profile directly
    # my ($x, $y, $xy) = @_;
    mutual_info = 0;
    for i in range(0, 21):
        for j in range(0, 21):
            if (xy[i][j] > 0):
                # comment:$xy->[$i][$j] > 0  ==> $x->[$i] > 0 && $y->[$j] > 0
                # print $xy->[$i][$j], ", ", $x->[$i],",",  $y->[$j], "\n";
                mutual_info += (xy[i][j] * math.log(xy[i][j] / (x[i] * y[j])))
    return mutual_info


# amino acid type info
def aatype(aa1, aa2):  # 10 type index: 0,1,2,3,4,5,6,7,8,9

    # my ($aa1, $aa2) = @_;
    nonpolar = "GAVLIPMFW"
    polar = "STNQCY"
    acidic = "DE"
    basic = "KRH"
    type1 = 0
    type2 = 0
    if (polar.find(aa1) >= 0):
        type1 = 1
    if (acidic.find(aa1) >= 0):
        type1 = 2
    if (basic.find(aa1) >= 0):
        type1 = 3
    if (polar.find(aa2) >= 0):
        type2 = 1
    if (acidic.find(aa2) >= 0):
        type2 = 2
    if (basic.find(aa2) >= 0):
        type2 = 3
    type3 = 0
    if (type1 > type2):
        tmp = type1
        type1 = type2
        type2 = tmp
    if (type1 == 0):
        type3 = type2
    elif (type1 == 1):
        type3 = 3 + type2
    elif (type1 == 2):
        type3 = 5 + type2
    else:
        type3 = 9
    return type3


def amino_type(aa1):  # four type: 0,1,2,3

    # my $aa1 = $_[0];
    nonpolar = "GAVLIPMFW"
    polar = "STNQCY"
    acidic = "DE"
    basic = "KRH"
    type1 = 0
    if (polar.find(aa1) >= 0):
        type1 = 1
    if (acidic.find(aa1) >= 0):
        type1 = 2
    if (basic.find(aa1) >= 0):
        type1 = 3
    return type1

# return value
# 1
