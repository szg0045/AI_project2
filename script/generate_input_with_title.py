#!/usr/bin/perl -w
import argparse
import os
import subprocess
import string
# from sklearn import *
import numpy
from scipy import spatial
from scipy.stats.stats import pearsonr
import math
from potential import aatype, mutual, getLevittCP, getJerniganCP, getBraunCP, entropy
from potential import round as round_pl

'''
#Generate inputs for one sequence
#input: script path, sequence, ss, sa, corr, alignment file, separation
#Output: SVM-light type output
'''


def main():
    '''
	"need 5 parameters: package path, input seq file, alignment file, separation, dist threshold.\n"
	'''

    # if (@ARGV != 5)
    # {
    #	die "need 5 parameters: package path, input seq file, alignment file, separation, dist threshold.\n";
    # }

    args = handle_args()
    if args['package_path'] is None:
        raise ValueError("The package_path directory is required.")
    if args['seq_file'] is None:
        raise ValueError("The input seq file must be provided.")
    if args['ali_file'] is None:
        raise ValueError("The alignment file must be provided.")
    if args['separation'] is None:
        raise ValueError("The separation file must be provided.")
    if args['threshold'] is None:
        raise ValueError("The dist threshold file must be provided.")
    # if args['input'] is None:
    #	raise ValueError("The input file must be provided.")
    # if args['output'] is None:
    #	raise ValueError("The output file must be provided.")
    # $package_path = shift @ARGV;
    package_path = str(args['package_path'])
    # $package_path='/home/project/conD/apps/SVMcon/svmcon1.0/script/'
    if (package_path[(len(package_path) - 1):] != "/"):
        package_path += "/"
    pot_address = package_path + "potential.py"
    # import pot_address
    # require "${package_path}potential.pl";
    seq_file = str(args['seq_file'])
    # $seq_file = shift @ARGV;
    # open(SEQ, "$seq_file") || die "can't read sequence file: $seq_file\n";
    # $seq_arg = <SEQ>; #separated by space
    # chomp $seq_arg;
    SEQ = open(seq_file)
    seq_arg = SEQ.readline()
    seq_arg = seq_arg[:-1]

    # $ss_arg = <SEQ>; #separated by space
    # chomp $ss_arg;
    ss_arg = SEQ.readline()  # separated by space
    ss_arg = ss_arg[:-1]
    # $sa_arg = <SEQ>; #separated by space (10-buried, 50 exposed)
    # chomp $sa_arg;
    sa_arg = SEQ.readline()  # separated by space (10-buried, 50 exposed)
    sa_arg = sa_arg[:-1]
    # $coor_arg = <SEQ>; #coordinates
    # chomp $coor_arg;
    coor_arg = SEQ.readline()  # coordinates
    # coor_arg = coor_arg[:-1]

    # close SEQ;

    # $ali_file = shift @ARGV;
    ali_file = str(args['ali_file'])
    # $separation = shift @ARGV; #Only predict contacts with separation equal or greater than 6.
    separation = str(args['separation'])
    # if ($separation < 6) { die "sequence separation is at least 6.\n"; };
    if (separation < 6):
        raise ValueError("sequence separation is at least 6.\n")
    # $threshold = shift @ARGV;  #contact threashold
    threshold = str(args['threshold'])

    # $win_size = 9; #size of window around aa of interest
    win_size = 9;  # size of window around aa of interest
    # $seg_size = 5; #size of window of central segment
    seg_size = 5;  # size of window of central segment

    # print seq_arg, len(seq_arg)
    # print ss_arg, len(ss_arg)
    # print sa_arg, len(sa_arg)
    # print coor_arg, len(coor_arg)


    # @seq = split(/\s+/, $seq_arg);
    seq = seq_arg.split()
    # @ss = split(/\s+/, $ss_arg);
    ss = ss_arg.split()
    # @sa = split(/\s+/, $sa_arg);
    sa = sa_arg.split()
    # @coor = split(/\s+/, $coor_arg);
    coor = coor_arg.split()

    # $length = @seq;
    length = len(seq)

    # if (@seq != @ss || @seq != @sa || @coor!=3*@seq)
    # {
    #	die "length is not consistent, $ali_file\n";
    # }
    if (len(seq) != len(ss) or len(seq) != len(sa) or len(coor) != 3 * len(seq)):
	import pdb; pdb.set_trace()
        raise ValueError("length is not consistent, $ali_file\n")
    # open(ALI, "$ali_file") || die "can't read $ali_file\n";
    # @ali = <ALI>;
    # close ALI;

    ALI = open(ali_file)
    ali = ALI.readlines()

    # generate profile
    # In future the MSA must be filtered or weighted by identity
    # $ali_num = shift @ali;
    ali_num = int(ali[0])
    # chomp $ali_num;
    ali = ali[1:]
    for i in range(0, len(ali)):
        ali_aa = ali[i]
        # if ali_aa[0] == '.':
        # 	ali_aa = ali_aa[1:]
        if ali_aa[-1] == '\n':
            ali_aa = ali_aa[:-1]
        ali[i] = ali_aa
    # print "length of ", i, " ", len(ali_aa)

    # print "ali: ", ali

    # $join_seq = join("", @seq);
    join_seq = "".join(seq)
    # push @ali, "$join_seq\n";
    ali.append(join_seq + '\n')
    # $ali_num++;
    # ali_num = ali_num + 1
    # $aa_str = "ACDEFGHIKLMNPQRSTVWY"; #20 std aa
    aa_str = "ACDEFGHIKLMNPQRSTVWY"  # 20 std aa

    # profile size: 20 + one gap, if all zero, means not exist(for pos outside of window)
    # @profile = ();
    profile = [[None] * 21] * length
    '''
	for ($i = 0; $i < $length; $i++)
	{
		for ($j = 0; $j < 21; $j++)
		{
			$profile[$i][$j] = 0;
		}
	}
	for ($i = 0; $i < $length; $i++)
	{
		for($j = 0; $j < $ali_num; $j++)
		{
			$aa = substr($ali[$j], $i, 1);
			$aa = uc($aa);
			$idx = index($aa_str, $aa);
			if ($idx < 0) #gap case or unknonw
			{
				#treated as a gap
				$idx = 20;
			}
			$profile[$i][$idx] +=  (1 / $ali_num);
		}
	}
	'''
    # print profile
    for i in range(0, length):
        for j in range(0, 21):
            profile[i][j] = 0
    #import pdb; pdb.set_trace()
    for i in range(0, length):
        for j in range(0, ali_num):
            aa = ali[j][i]
            aa = string.upper(aa)
            # print "aa: ", aa
            if aa == '.' or aa == 'X':
                idx = -1
            else:
                try:
                    idx = aa_str.find(aa)
                except:
                    import pdb;
                    pdb.set_trace()
            if (idx < 0):  # gap case or unknonw
                # treated as a gap
                idx = 20;
            profile[i][idx] += (1.0 / ali_num)

    # print profile
    # print length
    # print type(separation)
    # print separation
    for i in range(0, length):
        for j in range(0, i):  # Only generate for lower triangle for saving space
            if abs(i - j) < int(separation):
                # next;
                continue
            # import pdb; pdb.set_trace()
            # pairwise information
            # @input_pair = ();
            input_pair = []
            # input window around target
            # @input_win = ();
            input_win = []
            # central segment information
            # @input_seg = ();
            input_seg = []
            # average of central segment
            # @input_seg_ave = ();
            input_seg_ave = []
            # average of protein information
            # @input_prot = ();
            input_prot = []
            # combination of all inputs
            # @total_input = ();
            total_input = []
            #####################Pairwise Information###################################
            # generate pairwise information(right now only for target pair)
            # Can be extended to all pairs in window and can be average of profile.
            # @prof1 = ();
            prof1 = []
            # @prof2 = ();
            prof2 = []
            for k in range(0, 21):
                # push @prof1, $profile[$i][$k];
                prof1.append(profile[i][k])
                # push @prof2, $profile[$j][$k];
                prof2.append(profile[j][k])
            # $vcos = &cosine(\@prof1, \@prof2); #this one call function cosine similarity
            vcos = 1 - spatial.distance.cosine(prof1, prof2)
            # $vcor = &correlation(\@prof1, \@prof2);  #calculates a Pearson correlation coefficient
            vcor = numpy.corrcoef(prof1, prof2)[0, 1]
            # @vaat = (); #amino acid type(10 different combination)
            vaat = []  # amino acid type(10 different combination)
            for k in range(0, 10):
                # push @vaat, 0;
                vaat.append(0)
            amino1 = seq[i]
            amino2 = seq[j]
            # idx=aa_str.index(aa)
            if (aa_str.find(amino1) > -1 and aa_str.find(amino2) > -1):
                # try:
                amino_type = aatype(amino1, amino2)
                vaat[amino_type] = 1
            # except:
            # 	import pdb; pdb.set_trace()
            # generate joint distribution mutual information
            # shouldn't be a 21 * 21 matrix, almost half is redudant
            # either make them symetic or train two triangles later
            joint_dist = [[0] * 21] * 21
            # for k in range(0,21):
            # 	for m in range(0,21):
            # 		joint_dist[k][m] = 0
            for k in range(0, ali_num):
                # $aa1 = substr($ali[$k], $i, 1);
                aa1 = ali[k][i]
                # $aa1 = uc($aa1);
                aa1 = string.upper(aa1)
                # $aa2 = substr($ali[$k], $j, 1);
                aa2 = ali[k][j]
                # $aa2 = uc($aa2);
                aa2 = string.upper(aa2)
                # $idx1 = index($aa_str, $aa1);
                idx1 = aa_str.find(aa1)
                # $idx2 = index($aa_str, $aa2);
                idx2 = aa_str.find(aa2)
                if (idx1 < 0):
                    idx1 = 20
                if (idx2 < 0):
                    idx2 = 20
                joint_dist[idx1][idx2] += (1 / ali_num)
            # $joint_dist[$idx1][$idx2] += ( 1 / (2 * $ali_num));
            # $joint_dist[$idx2][$idx1] += ( 1 / (2 * $ali_num));
            # print "joint_dist: ", joint_dist
            mutual_info = mutual(prof1, prof2, joint_dist)
            # total input size of pair-wise: 13
            # push @input_pair, $vcos, $vcor, $mutual_info, @vaat;
            input_pair.append(vcos)
            input_pair.append(vcor)
            input_pair.append(mutual_info)
            # input_pair.append(vaat)
            input_pair = input_pair + vaat
            # get the three differnt types of pariwise potentials
            ##################End of pairwise information###############################
            # push @input_pair, &getLevittCP($amino1, $amino2);
            input_pair.append(getLevittCP(amino1, amino2))
            # push @input_pair, &getJerniganCP($amino1, $amino2);
            input_pair.append(getJerniganCP(amino1, amino2))
            # push @input_pair, &getBraunCP($amino1, $amino2);
            input_pair.append(getBraunCP(amino1, amino2))
            ##################Generate Window Information###############################
            # for each position in window(27): 21 for profile, 3SS, 2AA, 1 entropy
            # total input window size: 486
            # for ($k = 0; $k < 2 * $win_size * 27; $k++)
            # {
            #	$input_win[$k] = 0;
            # }
            input_win = [0] * 2 * win_size * 27
            # for k in range (0, 2*win_size*27):
            # 	input_win[k] = 0
            # $half_win = int($win_size / 2);
            half_win = int(win_size / 2)
            # create first window
            '''
			for ($k = $i - $half_win; $k <= $i + $half_win; $k++)
			{
				if ($k < 0 || $k >= $length) { next; };
				#copy profile
				$start = ($k - $i + $half_win ) * 27;
				for ($m = 0; $m < 21; $m++)
				{
					$input_win[$start + $m] = $profile[$k][$m];
				}
				#SS information
				$start += 21;
				$sec = $ss[$k];
				if ($sec eq "H") { $input_win[$start] = 1;}
				elsif ($sec eq "E") {$input_win[$start + 1] = 1; }
				elsif ($sec eq "C") {$input_win[$start + 2] = 1; };
				$start += 3;
				$sov = $sa[$k];
				if ($sov < 25) { $input_win[$start] = 1;}
				else { $input_win[$start + 1] = 1; };
				$start += 2;
				#Entropy information
				@prof1 = ();
				for ($m = 0; $m < 21; $m++)
				{
					push @prof1, $profile[$k][$m];
				}
				$input_win[$start] = &entropy(\@prof1);
			}
			#information for second window
			for ($k = $j - $half_win; $k <= $j + $half_win; $k++)
			{
				if ($k < 0 || $k >= $length) { next; };
				#copy profile
				$start = ($k - $j + $half_win + $win_size ) * 27;
				for ($m = 0; $m < 21; $m++)
				{
					$input_win[$start + $m] = $profile[$k][$m];
				}
				#SS information
				$start += 21;
				$sec = $ss[$k];
				if ($sec eq "H") { $input_win[$start] = 1;}
				elsif ($sec eq "E") {$input_win[$start + 1] = 1; }
				elsif ($sec eq "C") {$input_win[$start + 2] = 1; };
				$start += 3;
				$sov = $sa[$k];
				if ($sov < 25) { $input_win[$start] = 1;}
				else { $input_win[$start + 1] = 1; };
				$start += 2;
				#Entropy information
				@prof1 = ();
				for ($m = 0; $m < 21; $m++)
				{
					push @prof1, $profile[$k][$m];
				}
				$input_win[$start] = &entropy(\@prof1);
			}'''
            # create first window
            for k in range(i - half_win, i + half_win + 1):
                if (k < 0 or k >= length):
                    continue
                # copy profile
                start = (k - i + half_win) * 27
                for m in range(0, 21):
                    input_win[start + m] = profile[k][m]
                # SS information
                start += 21
                sec = ss[k]
                if (sec == "H"):
                    input_win[start] = 1
                elif (sec == "E"):
                    input_win[start + 1] = 1
                elif (sec == "C"):
                    input_win[start + 2] = 1
                start += 3
                sov = sa[k]
                if (sov < 25):
                    input_win[start] = 1
                else:
                    input_win[start + 1] = 1
                start += 2
                # Entropy information
                prof1 = []
                for m in range(0, 21):
                    # push @prof1, $profile[$k][$m];
                    prof1.append(profile[k][m])
                input_win[start] = entropy(prof1)
            # information for second window
            for k in range(j - half_win, j + half_win + 1):
                if (k < 0 or k >= length):
                    continue
                # copy profile
                start = (k - j + half_win + win_size) * 27
                for m in range(0, 21):
                    input_win[start + m] = profile[k][m]
                # SS information
                start += 21
                sec = ss[k]
                if (sec == "H"):
                    input_win[start] = 1
                elif (sec == "E"):
                    input_win[start + 1] = 1
                elif (sec == "C"):
                    input_win[start + 2] = 1
                start += 3
                sov = sa[k]
                if (sov < 25):
                    input_win[start] = 1
                else:
                    input_win[start + 1] = 1
                start += 2
                # Entropy information
                # @prof1 = ();
                prof1 = []
                for m in range(0, 21):
                    prof1.append(profile[k][m])
                input_win[start] = entropy(prof1)

            ##################End of Window Information#######################

            ##################Generate Segment Window Information#############
            # Total size of segment window is: 135
            input_seg = [0] * seg_size * 27
            # for k in range (0,seg_size * 27):
            # 	input_seg[k] = 0
            half_win = int(seg_size / 2)
            middle = int((i + j) / 2)
            # for ($k = $middle - $half_win; $k <= $middle + $half_win; $k++)
            for k in range(middle - half_win, middle + half_win):
                if (k < 0 or k >= length):
                    raise ValueError("should all be valid in segment")
                # copy profile
                start = (k - middle + half_win) * 27
                # for ($m = 0; $m < 21; $m++)
                for m in range(0, 21):
                    input_seg[start + m] = profile[k][m]
                # SS information
                start += 21
                sec = ss[k]
                if (sec == "H"):
                    input_seg[start] = 1
                elif (sec == "E"):
                    input_seg[start + 1] = 1
                elif (sec == "C"):
                    input_seg[start + 2] = 1
                start += 3
                sov = sa[k]
                if (sov < 25):
                    input_seg[start] = 1
                else:
                    input_seg[start + 1] = 1
                start += 2
                # Entropy information
                prof1 = []
                # for ($m = 0; $m < 21; $m++)
                for m in range(0, 21):
                    # push @prof1, $profile[$k][$m];
                    prof1.append(profile[k][m])
                input_seg[start] = entropy(prof1)

            # generate segment average information (Composition of SS, SA) and length information(in range or a number?)
            # is it composition of MSA or just target sequence???? (use profile instead of sequence)
            # 21 for composition of AA, 3 for SS, 11 for length (6,7,8,9,10-14,15-19,20-24,25-29,30-39,40-49,>49)
            # at this moment, use length (strictly separation:length+1) as
            # input (one input), so input size: 21 + 3 + 1 = 15.
            # Later, we might add composition of solvent accessibility

            # Total size of composition of segment: 42 (21 AA, 3 SS and 16 type, and 2 SA)
            input_seg_ave = [0] * 42
            # for k in range(0,42):
            # 	input_seg_ave[k] = 0
            seg_length = i - j - 1

            # WARNING: Here: We assume we work on Lower Triangle. When we work on Upper triangle, this
            # need to be changed.
            comp_bur = 0
            for k in range(j + 1, i - 1):  # ($k = $j + 1; $k <= $i - 1; $k++)
                # sum the profile
                for m in range(0, 21):  # ($m = 0; $m < 21; $m++)
                    input_seg_ave[m] += profile[k][m] / seg_length

                # sum the SS
                if (ss[k] == "H"):
                    input_seg_ave[21] += 1 / seg_length
                elif (seq[k] == "E"):
                    input_seg_ave[22] += 1 / seg_length
                elif (seq[k] == "C"):
                    input_seg_ave[23] += 1 / seg_length
                if (sa[k] < 25):
                    comp_bur += 1 / seg_length
            # set the separation between i and j ( = seg_length + 1)
            # more to do: add separation 10, 11, 12, 13. and check if casp separation include 12 or 24
            # usually they evaluate at: >=6, >=8, >=12, >=16, >=24, ....
            seg_length = seg_length + 1
            if (seg_length < 6):
                input_seg_ave[24] = 1
            elif (seg_length == 6):
                input_seg_ave[25] = 1
            elif (seg_length == 7):
                input_seg_ave[26] = 1
            elif (seg_length == 8):
                input_seg_ave[27] = 1
            elif (seg_length == 9):
                input_seg_ave[28] = 1
            elif (seg_length == 10):
                input_seg_ave[29] = 1  # added
            elif (seg_length == 11):
                input_seg_ave[30] = 1  # added
            elif (seg_length == 12):
                input_seg_ave[31] = 1  # added
            elif (seg_length == 13):
                input_seg_ave[32] = 1  # added
            elif (seg_length == 14):
                input_seg_ave[33] = 1
            elif (seg_length < 19):
                input_seg_ave[34] = 1
            elif (seg_length < 24):
                input_seg_ave[35] = 1
            elif (seg_length <= 29):
                input_seg_ave[36] = 1
            elif (seg_length <= 39):
                input_seg_ave[37] = 1
            elif (seg_length <= 49):
                input_seg_ave[38] = 1
            else:
                input_seg_ave[39] = 1
            input_seg_ave[40] = comp_bur
            input_seg_ave[41] = 1 - comp_bur
            ##################End of Segment Window Information###############

            ##################Generate Protein Composition Information############################
            # 21 for AA, 3 for SS, 4 for length, 2 for SA
            # Total size of protein information is: 30
            input_prot = [0] * 30
            # for k in range (0,30): #($k = 0; $k < 30; $k++)
            # 	input_prot[k] = 0
            for k in range(0, length):  # ($k = 0; $k < $length; $k++)
                # sum the profile
                for m in range(0, 21):  # ($m = 0; $m < 21; $m++)
                    input_seg_ave[m] += profile[k][m] / length
                # sum the SS
                if (ss[k] == "H"):
                    input_prot[21] += 1 / length
                elif (ss[k] == "E"):
                    input_prot[22] += 1 / length
                elif (ss[k] == "C"):
                    input_prot[23] += 1 / length
            # set the separation between i and j ( = seg_length + 1)
            ##################more do do: NEED TO CHANGE ACCODING TO TRAIN DATAET: better use: 50, 100, 150, >=200 above#########
            if (length <= 50):
                input_prot[24] = 1
            elif (length <= 100):
                input_prot[25] = 1
            elif (length <= 150):
                input_prot[26] = 1
            # elsif ($length <= 200) { $input_prot[27] = 1; }
            else:
                input_prot[27] = 1

            # composition of solvent acc

            for k in range(0, length):  # ($k = 0; $k < $length; $k++)
                if (sa[k] < 25):
                    input_prot[28] += 1 / length
                else:
                    input_prot[29] += 1 / length

            ##################End of Protien Composition Information##############################


            # round the numbers into at most 2 decimal
            # $total_size = 16 + 18 * 27[=486] + 5 * 27[=135] + 42 + 30 = 709;
            total_size = 709
            # print "pair: ", join(",", @input_pair), "\n";
            # print "win: ", join(",", @input_win), "\n";
            # $tmp_size = @input_seg;
            # print "segment size: $tmp_size\n";
            # print "seg: ", join(",", @input_seg) , "\n";
            # print "seg_ave: ", join(",", @input_seg_ave), "\n";
            # print "prot: ", join(",", @input_prot), "\n";
            #####push @total_input, @input_pair, @input_win, @input_seg, @input_seg_ave, @input_prot;
            # print "input pair: ", input_pair, len(input_pair)
            # print "input win: ", input_win, len(input_win)
            # print "input seg: ", input_seg, len(input_seg)
            # print "input seg ave: ", input_seg_ave, len(input_seg_ave)
            # print "input prot: ", input_prot, len(input_prot)
            total_input = input_pair + input_win + input_seg + input_seg_ave + input_prot
            # total_input = input_pair
            # total_input.append(input_win)
            # total_input.append(input_seg)
            # total_input.append(input_seg_ave)
            # total_input.append(input_prot)
            # print total_input, len(total_input)
            # $input_size = @total_input;
            # print "input size: $input_size\n";
            for k in range(0, total_size):  # (k= 0; $k < $total_size; $k++)
                # print "$total_input[$k]\n";
                total_input[k] = round_pl(total_input[k])

            # Generate label
            # $x1=$y1=$z1=$x2=$y2=$z2=0;
            x1 = y1 = z1 = x2 = y2 = z2 = 0
            x1 = float(coor[3 * i])
            y1 = float(coor[3 * i + 1])
            z1 = float(coor[3 * i + 2])
            x2 = float(coor[3 * j])
            y2 = float(coor[3 * j + 1])
            z2 = float(coor[3 * j + 2])
            dist = math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
            label = "-1"
            if (dist < threshold):
                label = "+1"
            # pint out the title comments (indices for residue pairs)
            print "#", j + 1, " ", i + 1

            # print label and inputs
            row = str(label)
            for k in range(0, total_size):  # ($k = 0; $k < $total_size; $k++)
                # if (total_input[k] != 0):
                row = row + " " + str(k + 1) + ":" + str(total_input[k])
            print row
        # print $label, " ", join(",", @total_input), "\n";
        # print "PRESS any key to continue....\n";

def handle_args():
    """
	Parses the command line arguments
	:return: args
	"""
    parser = argparse.ArgumentParser()
    # need 5 parameters:  package path, input seq file, alignment file, separation, dist threshold
    parser.add_argument("package_path", help="The package_path directory is required.")
    parser.add_argument("seq_file", help="The input seq file must be provided.")
    parser.add_argument("ali_file", help="The alignment file must be provided")
    parser.add_argument("separation", help="The separation file must be provided.")
    parser.add_argument("threshold", help="The dist threshold file must be provided.")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    main()
