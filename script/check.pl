#! /usr/bin/perl -w


$svm_predictor = shift @ARGV;
$input_file = shift @ARGV;
$svm_model = shift @ARGV;
$output_file = shift @ARGV;
$svm_predictor = '/Users/saurabhgupta/Documents/script/server/svm_classify';
$svm_model = '/Users/saurabhgupta/Documents/script/model/model.g3';
if (! -f $svm_predictor)
{
	die "can't find svm classifier.\n";
}
if (! -f $svm_model)
{
	die "can't find the model definition file.\n";
}
if (! -f $input_file)
{
	die "can't find the fasta file.\n";
}

print "classify data points using SVM(perl)...\n";
#make svm predictions
system("$svm_predictor $input_file $svm_model $output_file");

print "Generated output_file"
