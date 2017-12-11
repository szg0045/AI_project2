import os


def main():
    script_dir = os.getcwd()
    print script_dir
    tmp_file = 'output_file1aqta.tmp'
    input_file = '/1aqta.seq'
    output_file = 'output_file1aqta'
    os.system("python {0} '{1}' {2} '{3}' {4}".format(
        'predict_map.py',
        script_dir,
        tmp_file,
        script_dir + input_file,
        output_file,
    )
    )


if __name__ == "__main__":
    main()
