import sys, itertools, re

sciteFile = open(sys.argv[1])
sciphiFile = open(sys.argv[2], "w+")

seq = []
for l in sciteFile:
        x = l.strip().split()
        for i, v in enumerate(x):
                while i >= len(seq):
                        seq.append([])
                seq[i].append(int(v))

def sciteChar2Phiscs(c):
	if c == 0:
		return '0'
	elif c == 1:
		return '1'
	elif c == 3:
		return '?'

if len(seq) > 0:
	print('cellID/mutID\t{}'.format('\t'.join(['mut'+str(i) for i in range(len(seq[0])) ])), file=sciphiFile)
	for i, s in enumerate(seq):
		print('cell{}\t{}'.format(i+1, '\t'.join([sciteChar2Phiscs(c) for c in s])), file=sciphiFile)
