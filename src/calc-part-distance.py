import sys
import scipy.optimize as opt
import numpy as np

def loadPart(fName):
	f = open(fName)
	ret = []
	for l in f:
		ret.append(l.strip().split())
	return ret

def eqDistance(P1, P2, normalize):
	P1 = [p for p in P1 if len(p) > 1]
	P2 = [p for p in P2 if len(p) > 1]
	P1.sort()
	P2.sort()
	i = j = common = 0
	while i < len(P1) and j < len(P2):
		if P1[i] == P2[j]:
			i += 1
			j += 1
			common += 1
		elif P1[i] < P2[j]:
			i += 1
		else:
			j += 1
	normalFactor = 1
	if normalize != 0:
		if normalize == 1:
			normalFactor = len(P1)
		if normalize == 2:
			normalFactor = len(P2)
	return common / normalFactor

def matchDistance(P1, P2, normalize):
	P1 = [p for p in P1 if len(p) > 1]
	P2 = [p for p in P2 if len(p) > 1]
	cost = np.zeros((len(P1), len(P2)))
	for i in range(len(P1)):
		for j in range(len(P2)):
			v = len(set(P1[i]).symmetric_difference(P2[j]))
			v = max(v, allCellCount - v)
			cost[i][j] = -v
	#print(cost, file=sys.stderr)
	row_ind, col_ind = opt.linear_sum_assignment(cost)
	#print("{}, {}: {}, {} : {} // {} ".format(len(P1), len(P2), len(row_ind), len(col_ind), row_ind, col_ind), file=sys.stderr)
	normalFactor = 1
	if normalize != 0:
		if normalize == 1:
			normalFactor = len(P1) * allCellCount
		if normalize == 2:
			normalFactor = len(P2) * allCellCount
		if normalFactor == 0: 
			normalFactor = 1
	return -sum([cost[row_ind[i], col_ind[i]] for i in range(min(len(P1), len(P2)))]) / normalFactor
	
argFirstListIndex = 1
distanceMethod = eqDistance
if sys.argv[argFirstListIndex] == '--match':
	distanceMethod = matchDistance
	argFirstListIndex += 1
	
normalize = 0
if sys.argv[argFirstListIndex] == '--normalize' or sys.argv[argFirstListIndex] == '--normalize-r':
	if sys.argv[argFirstListIndex] == '--normalize':
		normalize = 1
	if sys.argv[argFirstListIndex] == '--normalize-r':
		normalize = 2
	argFirstListIndex += 1
	

T, P = [], []
for i in range(argFirstListIndex, len(sys.argv), 2):
	fn = sys.argv[i]
	p = loadPart(fn)
	P.append(p)
	T.append(sys.argv[i+1])

allCells = set.union(*[set.union(*[set(pp) for pp in p]) for p in P])
allCellCount = len(allCells)

print("{}".format(' '.join([t for t in T])))
for i in range(len(T)):
	print(T[i], end=" ")
	for j in range(len(T)):
		#print(P[i])
		d = distanceMethod(P[i], P[j], normalize)
		print("{0}".format(d), end=" ")
	print()
