#!/bin/sh
if [ $# -ne 3 ]
then
	echo "need three parameters: data set, model, output file."
	exit 1
fi
/Users/saurabhgupta/Documents/script/server/svm_classify $1 $2 $3
