import sys, itertools, re

def distance(seq, mutationSet):
	# print(len(seq))
	# print(seq)
	# print(mutationSet)
	s2 = [0 for s in seq]
	for m in mutationSet:
		if m-1 < len(seq):
			s2[m-1] = 1
	d = 0
	for x, y in zip(seq, s2):
		if x != y and x != 3:
			d += 1
	return d

sciteFile = open(sys.argv[1])
seqFile = open(sys.argv[2])
treeFile = open(sys.argv[3], "w+")
cloneFile = open(sys.argv[4], "w+")

sequences = []
for l in seqFile:
	for num, val in enumerate(l.strip().split()):
		while num >= len(sequences):
			sequences.append([])
		sequences[num].append(int(val))


muts = []
edges = []
parent = {}
clones = {}

inGraph = False
for l in sciteFile:
	m = re.match("digraph g {", l)
	if m:
		inGraph = True
		continue
	if inGraph and re.match("}", l):
		inGraph = False
	if inGraph:
		m = re.search('"(\d+)" -> "(\d+)";', l)
		if m:
			v, u = int(m.group(1)), int(m.group(2))
			edges.append((v,u))
			muts.append(v)
			muts.append(u)
			parent[u] = v
		m = re.search('"(\d+)" \[label="(\d+)"\];', l)
		if m:
			#print("clone found {} {}".format(v, u))
			v, u = int(m.group(1)), int(m.group(2))
			clones[v] = [u]


allMutations = {}
mark = {}

def dfs(m):
	mark[m] = True
	allMutations[m] = []
	if m in clones:
		allMutations[m] = [clones[m][0]]
	if m not in parent:
		return
	if parent[m] not in mark:
		dfs(parent[m])
	allMutations[m] = allMutations[parent[m]] + allMutations[m]

for m in muts:
	if m not in mark:
		dfs(m)


muts = set(muts)
print(' '.join([str(c) for c in muts]), file = treeFile)
for (v, u) in edges:
	# print("{} -> {}".format(v, u))
	# Order of vertices in our tree is reverse of the order of vertices of SCITE
	d = len(set(allMutations[v]).symmetric_difference(allMutations[u]))
	print("{}->{} {}".format(u, v, d), file=treeFile)
	#print("{}->{} {}".format(u, v, ','.join([str(c) for c in allMutations[v]])), file=treeFile)

for m in muts:
	if m not in clones:
		clones[m] = []

for cl, cells in clones.items():
	print("{0}".format(cl), end=" ", file=cloneFile)
	for c in cells:
		print("{0}".format(c), end = " ", file=cloneFile)
	print(file=cloneFile)



