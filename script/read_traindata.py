import numpy
import pdb
import re
def read_traindata(filename):#you should provide the filename
	#filename = "train_data.txt"
	#train_data=numpy.zeros((267793, 709))
	train_data=[]
	train_label=[]#numpy.zeros((267793, 1))

	#for word in f.read().split():
	#    print(word)
	i=0
	with open(filename) as f:
	    for line in f:
		train_data.append([0 for k in range(0,709)])
		splitted_line=re.split('[ :]+', line.strip())
		#train_label[i]=float(splitted_line[0])#float(line.strip().split(" ")[0])
		train_label.append(float(splitted_line[0]))
		splitted_line.pop(0)	
		for k in range(0,len(splitted_line)):
				if (k%2==0):
	        			#train_data[i,(int)(splitted_line[k])-1]=(float)(splitted_line[k+1])
					train_data[i][(int)(splitted_line[k])-1]=(float)(splitted_line[k+1])
		i=i+1
		if (i%1000==0):
			print("read train data: i=",i)
	return train_data,train_label
	
#print (train_label)
#print (train_data[2])

#pdb.set_trace()

