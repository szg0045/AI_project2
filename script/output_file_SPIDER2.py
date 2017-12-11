#!/usr/bin/env python
from sys import *

# Read .spd3 file and extract aa and ss and store them row wise in output_file.tmp and put ' ' in between each entry.

f=open('/home/project/conD/data/1aqta.spd3','r')
g=open('/home/project/conD/data/output_file.tmp','w')
next(f)
length=0
for line in f:
	parts=line.split()
	g.write(parts[1]+' ')
	length+=1
g.write('\n')
f.close()
g.close()

	
f=open('/home/project/conD/data/1aqta.spd3','r')
g=open('/home/project/conD/data/output_file.tmp','a')
next(f)
for line in f:
	parts=line.split()
	g.write(parts[2]+' ')
g.write('\n')
f.close()
g.close()

# ASA converter
dict_ASA0 = dict(zip("ACDEFGHIKLMNPQRSTVWY",
					(115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
		185, 160, 145, 180, 225, 115, 140, 155, 255, 230)))
if __name__ == '__main__':
	#if len()<2:
	#	print >>stderr, 'usage: RUN *.spd3'
	#	exit(1)
	g=open('/home/project/conD/data/output_file.tmp','a')
	for x in open('/home/project/conD/data/1aqta.spd3','r'):
		if x[0] == '#': continue
		ss = x.split()
		rASA = float(ss[3]) / dict_ASA0[ss[1]] * 100
		if rASA<=25.0:
			rASA=50
		else:
			rASA=10
		g.write(str(rASA)+' ')
		#print ss[1], '%.2f' % rASA
	g.write('\n')
g.close()

##Place holder for coordinates to generate label

g=open('/home/project/conD/data/output_file.tmp','a')

for i in range(length):
	if i<length-1:
		g.write('0 0 0 ')
	else:
		g.write('0 0 0\n')

g.close()
