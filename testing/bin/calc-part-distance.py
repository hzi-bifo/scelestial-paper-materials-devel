import sys, os
import scipy.optimize as opt
import numpy as np
from scipy.special import comb

def loadPart(fName):
	if not os.path.exists(fName):
		return None
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

def rfDistance(P1, P2, normalize):
	eqDist = eqDistance(P1, P2, 0)
	dist = len(P1) - eqDist + len(P2) - eqDist
	normalFactor = 1
	if normalize != 0:
		if normalize == 1:
			normalFactor = len(P1)
		if normalize == 2:
			normalFactor = len(P2)
	return dist / normalFactor

def mismatchDistance(P1, P2, normalize):
	P1 = [p for p in P1 if len(p) > 1]
	P2 = [p for p in P2 if len(p) > 1]
	cost = np.zeros((len(P1), len(P2)))
	for i in range(len(P1)):
		for j in range(len(P2)):
			#v = len(set(P1[i]).symmetric_difference(P2[j]))
			#v = max(v, allCellCount - v) 
			#The best of two cases is min, Cases: A <-> B, A' <-> B' and A' <-> B, A <-> B'
			A, B = set(P1[i]), set(P2[j])
			v = len(A.intersection(B)) + allCellCount - len(A.union(B))
			u = len(set(A).difference(B)) + len(set(B).difference(A))
			cost[i][j] = min(v, u)
			#print('{}:{}/{}'.format(cost[i][j], v, u), end = ' ', file=sys.stderr)
		#print(file=sys.stderr)
	#print(cost, file=sys.stderr)
	row_ind, col_ind = opt.linear_sum_assignment(cost)
	#print("{}, {}, {}: {}, {} : {} // {} ".format(len(P1), len(P2), allCellCount, len(row_ind), len(col_ind), row_ind, col_ind), file=sys.stderr)
	normalFactor = 1
	if normalize != 0:
		if normalize == 1:
			normalFactor = len(P1) * allCellCount
		if normalize == 2:
			normalFactor = len(P2) * allCellCount
		if normalFactor == 0: 
			normalFactor = 1
	return sum([cost[row_ind[i], col_ind[i]] for i in range(min(len(P1), len(P2)))]) / normalFactor

def matchDistance(P1, P2, normalize):
	P1 = [p for p in P1 if len(p) > 1]
	P2 = [p for p in P2 if len(p) > 1]
	cost = np.zeros((len(P1), len(P2)))
	for i in range(len(P1)):
		for j in range(len(P2)):
			#v = len(set(P1[i]).symmetric_difference(P2[j]))
			#v = max(v, allCellCount - v) 
			#The best of two cases is min, Cases: A <-> B, A' <-> B' and A' <-> B, A <-> B'
			A, B = set(P1[i]), set(P2[j])
			v = len(A.intersection(B)) + allCellCount - len(A.union(B))
			u = len(set(A).difference(B)) + len(set(B).difference(A))
			cost[i][j] = -max(v, u)
			#print('{}:{}/{}'.format(cost[i][j], v, u), end = ' ', file=sys.stderr)
		#print(file=sys.stderr)
	#print(cost, file=sys.stderr)
	row_ind, col_ind = opt.linear_sum_assignment(cost)
	#print("{}, {}, {}: {}, {} : {} // {} ".format(len(P1), len(P2), allCellCount, len(row_ind), len(col_ind), row_ind, col_ind), file=sys.stderr)
	normalFactor = 1
	if normalize != 0:
		if normalize == 1:
			normalFactor = len(P1) * allCellCount
		if normalize == 2:
			normalFactor = len(P2) * allCellCount
		if normalFactor == 0: 
			normalFactor = 1
	return -sum([cost[row_ind[i], col_ind[i]] for i in range(min(len(P1), len(P2)))]) / normalFactor

def subprobDistance(P1, P2, normalize):
	def partProb(allCellCount, a, b, iab):
		choose = comb
		r = 0
		for i in range(iab, a+1):
			r += choose(allCellCount-i, b-i) * choose(a, i)
		return r/choose(allCellCount, b)
	def completePart(P2):
		P2.extend([allCells.difference(set(p)) for p in P2 if len(p) > 1])
		return P2
#print("P1:{}, P2: {}".format(P1, P2))
	P1 = completePart(P1)
	P2 = completePart(P2)
	r = []
	for i in range(len(P1)):
		min_p, min_j, min_B = 100, -1, None
		for j in range(len(P2)):
			#v = len(set(P1[i]).symmetric_difference(P2[j]))
			#v = max(v, allCellCount - v) 
			#The best of two cases is min, Cases: A <-> B, A' <-> B' and A' <-> B, A <-> B'
			A, B = set(P1[i]), set(P2[j])
			p = partProb(allCellCount, len(A), len(B), len(A.intersection(B))) 
			if p < min_p:
				min_p, min_j, min_B = p, j, B
		r.append((A, min_B, min_p))
	return r


argFirstListIndex = 1
distanceMethod = eqDistance
#the default case
if sys.argv[argFirstListIndex] == '--eq':
	distanceMethod = eqDistance
	argFirstListIndex += 1
	
if sys.argv[argFirstListIndex] == '--rf':
	distanceMethod = rfDistance
	argFirstListIndex += 1
	
if sys.argv[argFirstListIndex] == '--match':
	distanceMethod = matchDistance
	argFirstListIndex += 1
	
if sys.argv[argFirstListIndex] == '--mismatch':
	distanceMethod = mismatchDistance
	argFirstListIndex += 1

if sys.argv[argFirstListIndex] == '--tree-prob':
	distanceMethod = subprobDistance
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

#allCells = set.union(*[set.union(*[set(pp) for pp in p]) for p in P])
allCells = set([j for pp in P if pp is not None for i in pp for j in i])
#print(allCells)
allCellCount = len(allCells)

print("{}".format(' '.join([t for t in T])))
for i in range(len(T)):
	print(T[i], end=" ")
	for j in range(len(T)):
		#print(P[i])
		#print('{},{} / {}, {}'.format(T[i], T[j], len(P[i]), len(P[j])), file=sys.stderr)
		if P[i] is not None and P[j] is not None:
			d = "{0}".format(distanceMethod(P[i], P[j], normalize))
		else:
			d = 'NA'
		print(d, end=" ")
	print()
