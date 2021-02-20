import sys, itertools, re

sciphiFile = open(sys.argv[1])
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
for l in sciphiFile:
	m = re.search("(\\d+)->(\\d+) ;", l)
	if m:
		v, u = int(m.group(1)), int(m.group(2))
		edges.append((v,u))
		muts.append(v)
		muts.append(u)
		parent[u] = v
	m = re.search("(\\d+)[[]shape=box,style=filled, fillcolor=white,label=\"([^\"]*)\"];", l)
	if m:
		cl, cellNames = int(m.group(1)), m.group(2)
		cells = [int(re.search("PAT(\\d+)", c).group(1)) for c in cellNames.split(",")]
		clones[cl] = cells

muts = set(muts)
for c in muts:
	print(c, file=treeFile, end = " ")
print(file = treeFile)

for c in muts:
	if c not in clones:
		clones[c] = []

for (v, u) in edges:
	# print("{} -> {}".format(v, u))
	# Order of vertices in our tree is reverse of the order of vertices of SCITE
	print("{}->{} {}".format(v, u, 1), file=treeFile)

	
for cl, cells in clones.items():
	print("{0}".format(cl), end=" ", file=cloneFile)
	for c in cells:
		print("{0}".format(c), end = " ", file=cloneFile)
	print(file=cloneFile)

