import numpy
import pdb
import re
from read_traindata import read_traindata
from read_testdata import read_testdata
from sklearn import svm
import pickle

filename1="train_data.txt"
filename2="output_file.svm"
train_data,train_label=read_traindata(filename1)
test_data,test_label=read_testdata(filename2)

#print "train data label: ",train_label
#print "test data label: ",test_label
#print "train data [0]: ",train_data[0]
#print "test data [0]: ",test_data[0]
clf = svm.SVR()
clf.fit(train_data,train_label) 
test_result=clf.predict(test_data)
print sum ( array ( test_result ) != array ( test_label ) ) *1.0/ len ( test_label )
#print "test result [0]: ",test_result[0:20]


