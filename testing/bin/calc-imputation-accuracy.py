import sys, os
from sys import argv

def loadMatrix(fn):
	r = []
	with open(fn) as f:
		for l in f:
			r.append([int(v) for v in l.strip().split(' ')])
	assert(all([len(rr) == len(r[0]) for rr in r]))
	return r

def dist(m1, m2):
	r = 0
	for m1r, m2r in zip(m1, m2):
		for a, b in zip(m1r, m2r):
			if a != b: 
				r += 1
				#print('D {} {}'.format(a, b), file=sys.stderr)
			else:
				#print('E {} {}'.format(a, b), file=sys.stderr)
				pass
	return r


argc = 1

fileInput, fileTrue = argv[argc], argv[argc+1]
argc += 2
matrixInput = loadMatrix(fileInput)
matrixTrue = loadMatrix(fileTrue)


for i in range(argc, len(sys.argv), 2):
	fn = sys.argv[i]
	matrix = loadMatrix(fn)
	d = dist(matrix, matrixTrue)
	dmax = dist(matrixInput, matrixTrue)
	print('name:{} {} {}'.format(sys.argv[i+1], d, dmax))
	
	
