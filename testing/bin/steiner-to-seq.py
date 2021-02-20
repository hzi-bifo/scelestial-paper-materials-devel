import sys, itertools, re

seqSciteFile = open(sys.argv[1])
steinerFile = open(sys.argv[2])


## Read the imput output file 
n = int(steinerFile.readline().strip())
treeNodes = []
cells = []
maxAc = -1
seq = []
for i in range(n):
	l = steinerFile.readline().strip()
	x = l.split()
	treeNodes.append(x[0])
	if x[1] == '1':
		cells.append(x[0])
	ac = sum([1 for s in x[2] if s == 'A'])
	# print("  mx: ac={} len={} x[2]={} i={}".format(ac, len(x[2]), x[2], i))
	if ac >= maxAc:
		maxAc = ac
		maxAcLen = len(x[2])
		treeRootSeqIdx = int(x[0])
	
	# if it is an actual (and not internal) node
	if x[1] == '1':
		if sum([1 for s in x[3] if s not in set(['A', 'C'])]) > 0:
			print('imputed sequence contains invalid char {}'.format(x[3]), file=sys.stderr)
			raise Exception('Invalid imputation')
		seq.append(x[3])

# print("max: ac:{} ac-len:{} idx:{} ".format(maxAc, maxAcLen, treeRootSeqIdx))

if len(seq) == 0: exit()
for j in range(len(seq[0])):
	print(' '.join(map(lambda v: str(v), [0 if s[j] == 'A' else 1 for s in seq])))

