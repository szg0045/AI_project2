import os


filename = '/Users/saurabhgupta/Documents/COMP 6600 - Artificial Intelligence/Final Project/code/script/model/svmcon_train.dataset'
script_dir = '/Users/saurabhgupta/Documents/COMP 6600 - Artificial Intelligence/Final Project/code/script'
align_file = "1aqta.ssaalign"
with open(filename) as f:
    first_line = f.readline()
    instances = int(first_line.split()[0])
    for i in range(instances):
        one = f.readline()    #discard 1st
        two = f.readline()    #discard 2nd
        seq = f.readline()
        seq = seq[:-1]
        ss = f.readline()  # separated by space
        ss = ss[:-1]
        sa = f.readline()  # separated by space (10-buried, 50 exposed)
        sa = sa[:-1]
        seven = f.readline()    #discard 7rd
        coor = f.readline()  # coordinates
        eight = f.readline()    #discard 8th
        with open('temp.txt', 'w') as f2:
            f2.write(seq)
            f2.write(ss)
            f2.write(sa)
            f2.write(coor)
        import pdb; pdb.set_trace()
        os.system("python '{0}' '{1}' {2} {3} {4} {5} >> {6}".format(
            script_dir + '/generate_input_with_title.py',
            script_dir,
            'temp.txt',
            align_file,
            '6',
            '8',
            'training_dataset.svm')
        )
        os.remove('temp.txt')
