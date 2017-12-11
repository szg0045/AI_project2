import argparse
import os
import subprocess
import math


def main():
    '''
	need 4 parameters: script dir , SPIDER2 predictor , input file in FASTA format, output file.
	'''
    args = handle_args()
    if args['script'] is None:
        raise ValueError("The Script directory is required.")
    if args['ssa'] is None:
        raise ValueError("The SPIDER2 predictor must be provided.")
    # if args['svm'] is None:
    #     raise ValueError("The SVM classifier must be provided.")
    # if args['model'] is None:
    #     raise ValueError("SVM trained model must be provided.")
    if args['input'] is None:
        raise ValueError("The input file must be provided.")
    if args['output'] is None:
        raise ValueError("The output file must be provided.")

    script_dir = str(args['script'])
    ssa_predictor = str(args['ssa'])
    svm_predictor = script_dir + "/server/svm_classify"
    svm_model = script_dir + "/model/model.g3"
    fasta_file = str(args['input'])
    output_file = str(args['output'])

    if not os.path.isdir(script_dir):
        raise ValueError("can't find script dir.")
    if not os.path.isfile(ssa_predictor):
        raise ValueError("can't find the SPIDER2 predictor.")
    if not os.path.isfile(svm_predictor):
        raise ValueError("can't find svm classifier.")
    if not os.path.isfile(svm_model):
        raise ValueError("can't find the model definition file.")
    if not os.path.isfile(fasta_file):
	print fasta_file
        raise ValueError("can't find the fasta file.")

    target_name = None
    try:
        with open(fasta_file) as f:
            target_name = f.read()
    except Exception as e:
        print e.args
        raise ValueError("can't open the fasta file.")

    target_name = target_name[1:]
    print 'target_name: ', target_name

    # generate alignment, predict ss and sa
    pos = fasta_file.rfind('/')
    print 'pos: ', pos
    if pos < 0:
        ssa_file = fasta_file
    else:
        ssa_file = fasta_file[pos + 1:]

    # nothing is achieved with ssa_file
    # validate regex ssa_file find this out $ssa_file =~ s/\.[^.]+$//;

    # ssa_file = os.path.splitext(ssa_file)[0] + ".ssa"
    print ssa_file
    align_file = "1aqta.ssaalign"

    print "predict secondary structure and solvent accessibility..."
    # # find this out `$ssa_predictor $fasta_file`;
    # # notice: two files are generated from ssa predictor: one is ssa output, one is alignment file.
    #
    # seq = None
    # ss = None
    # sa = None
    # length = None
    with open(ssa_file) as infile:
        seq = infile.readline()
        # chomp
        ss = infile.readline()
        # chomp
        sa = infile.readline()
    #     # chomp
    #
    length = len(seq)
    # if length != len(ss) or length != len(sa):
    #     raise ValueError("sequence length doesn't match.")
    #
    # with open(output_file, 'w') as outfile:
    #     for i in range(0, length):
    #         outfile.write(seq[i])
    #         if i < length - 1:
    #             outfile.write(" ")
    #         else:
    #             outfile.write("\n")
    #
    #     for i in range(0, length):
    #         outfile.write(ss[i])
    #         if i < length - 1:
    #             outfile.write(" ")
    #         else:
    #             outfile.write("\n")
    #
    #     for i in range(0, length):
    #         if sa[i] == 'e':
    #             outfile.write(str(50))
    #         else:
    #             outfile.write(str(10))
    #         if i < length - 1:
    #             outfile.write(" ")
    #         else:
    #             outfile.write("\n")
    #
    #     for i in range(0, length):
    #         outfile.write("0 0 0")
    #         if i < length - 1:
    #             outfile.write(" ")
    #         else:
    #             outfile.write("\n")

    # generate svm dataset for separation >=6
    print "generate SVM dataset..."
    if not os.path.isfile(output_file + '.svm'):
        os.system("python '{0}' '{1}' {2} {3} {4} {5} > {6}".format(
            script_dir + '/generate_input_with_title.py',
            script_dir,
            ssa_predictor,
            align_file,
            '6',
            '8',
            output_file + '.svm'
            )
        )
    # with open(output_file + '.svm') as f:
    #     subprocess.call([script_dir + "/generate_input_with_title.pl", script_dir, output_file + ".tmp",
    #                      align_file, "6", "8"], stdout=f, stderr=f)

    print "classify data points using SVM...\n"
    if not os.path.isfile(output_file + '.res'):
        print svm_predictor
        # make svm predictions
        os.system("./check.pl '{0}' {1} '{2}' {3}".format(
            svm_predictor,
            output_file + '.svm',
            svm_model,
            output_file + '.res',)
        )
    import pdb; pdb.set_trace()
    # subprocess.call([svm_predictor, output_file + ".svm", svm_model, output_file + ".res"],
    #                 stdout=f, stderr=f)

    scores = []
    with open(output_file + '.res') as resfile:
        scores = resfile.readlines()

    svm_set = []
    with open(output_file + '.svm') as svmfile:
        svm_set = svmfile.readlines()

    # generate CASP format file
    casp = output_file + ".casp"
    with open(casp, 'w') as caspfile:
        caspfile.write("PFRMAT RR\n")
        caspfile.write("TARGET ")
        caspfile.write(target_name)
        caspfile.write("AUTHOR SVMcon\n")
        caspfile.write("METHOD SVM contact map predictor (separation >= 6)\n")
        caspfile.write("MODEL  1\n")

        for i in range(1, length + 1):
            caspfile.write(seq[i - 1])
            if i % 50 == 0 or i == length:
                caspfile.write("\n")

        # check consistency
        if len(svm_set) != 2 * len(scores):
            raise ValueError("number of data points doesn't match with number of scores.")

        pairs = []
        for i in range(0, len(scores)):
            score = float(scores[i])
            pair = svm_set[i]
            if pair.endswith("\n"):
                pair = pair[:-len("\n")]
            id1, id2 = map(int, pair.split())

            if score > 0:
                sigmoid = 1 / (1 + math.exp(-1 * score))
                sigmoid *= 100000
                sigmoid = int(sigmoid)
                sigmoid /= 100000
                pairs.append({'id1': id1, 'id2': id2, 'score': sigmoid})

        # sort the pairs and out pairs here
        sorted_pairs = sorted(pairs, key=lambda k: k['id2'])
        out_pairs = sorted(sorted_pairs, key=lambda k: k['id1'])

        for i in range(0, len(out_pairs)):
            caspfile.write(out_pairs[i][0])
            caspfile.write(out_pairs[i][1])
            caspfile.write(" 0  8 ")
            caspfile.write(out_pairs[i][2])
            caspfile.write("\n")

        caspfile.write("END\n")

    # remove temporary files
    # os.remove(output_file + '.tmp')
    # os.remove(output_file + '.res')
    # os.remove(output_file + '.svm')
    # os.remove(ssa_file)
    # os.remove(align_file)


def handle_args():
    """
	Parses the command line arguments
	:return: args
	"""
    parser = argparse.ArgumentParser()
    parser.add_argument("script", help="The Script directory.")
    parser.add_argument("ssa", help="SPIDER2 predictor.")
    # parser.add_argument("svm", help="svm classifier.")
    # parser.add_argument("model", help="svm model (trained SVM model).")
    parser.add_argument("input", help="The input file in FASTA.")
    parser.add_argument("output", help="The output file.")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    main()
