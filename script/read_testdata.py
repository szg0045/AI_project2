import numpy as np
import pdb
import re

def read_testdata(filename):	
	#filename = "output_file.svm"
	#print filename
	#test_data=numpy.zeros((267793, 709))
	#test_label=numpy.zeros((267793, 1))
	test_data=[]
	test_label=[]

	#for word in f.read().split():
	#    print(word)
	i=0
	with open(filename) as f:
		for line in f:
			if not(re.match(r'^#', line)):
				test_data.append([0 for k in range(0,709)])
				splitted_line=re.split('[ :]+', line.strip())
				#train_label[i]=float(splitted_line[0])#float(line.strip().split(" ")[0])
				#test_label.append([float(splitted_line[0])])#float(line.strip().split(" ")[0])
				test_label.append(float(splitted_line[0]))
				splitted_line.pop(0)	
				for k in range(0,len(splitted_line)):
						if (k%2==0):
							#print(test_data)
	        					test_data[i][(int)(splitted_line[k])-1]=(float)(splitted_line[k+1])
				i=i+1
				if (i%1000==0):
					print("read test data: i=",i)
		#print (np.asarray(test_label))
		#print (np.asarray(test_data))
		#pdb.set_trace()	
		
			
	#return np.asarray(test_data),np.asarray(test_label)
	return test_data,test_label
#read_testdata("output_file.svm")
	
#print (train_label)
#print (train_data[2])

#pdb.set_trace()




