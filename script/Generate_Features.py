import argparse
import os
import subprocess
import os
import sys

def main():
    '''
	need 4 parameters: script dir , SPIDER2 predictor , input file in FASTA format, output file.
	'''
    #print "1_AI"
    args = handle_args()
    print "1_AI"
    if args['script'] is None:
        raise ValueError("The Script directory is required.")
    if args['ssa'] is None:
        raise ValueError("The SPIDER2 predictor must be provided.")
    if args['input'] is None:
        raise ValueError("The input file must be provided.")
    if args['output'] is None:
        raise ValueError("The output file must be provided.")

    script_dir = str(args['script'])
    ssa_predictor = str(args['ssa'])
    fasta_file = str(args['input'])
    output_file = str(args['output'])

    if not os.path.isdir(script_dir):
        raise ValueError("can't find script dir.")
    if not os.path.isfile(ssa_predictor):
        raise ValueError("can't find the SPIDER2 predictor.")
    if not os.path.isfile(fasta_file):
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

    ssa_file += ".spd3"
    print ssa_file

    print "predict secondary structure and solvent accessibility..."
    # find this out `$ssa_predictor $fasta_file`;
    # notice: two files are generated from ssa predictor: one is ssa output, one is alignment file.

    print "DeepCNF dataset..."
    # if not os.path.exists(output_file):
    os.system("python '{0}' '{1}' {2} {3} {4} {5} > {6}".format(
        script_dir + '/generate_input_with_title.py',
        script_dir,
        'output_file.tmp',
        '1aqta.ssaalign',
        '6',
        '8',
        output_file)
    )
    # with open(output_file) as f:
    #     #subprocess.call(["perl" +" " + script_dir + "/generate_input_with_title.pl",script_dir,"output_file.tmp",
    #      #                "1aqta.ssaalign", "6", "8"], stdout=f, stderr=f)
    #
    #     f.write(output)

def handle_args():
    """
	Parses the command line arguments
	:return: args
	"""
    parser = argparse.ArgumentParser()
    #print parser
    parser.add_argument("script", help="The Script directory.")
    parser.add_argument("ssa", help="SPIDER2 predictor.")
    parser.add_argument("input", help="The input file in FASTA.")
    parser.add_argument("output", help="The output file.")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    main()
