import sys, itertools, re

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


edges = []
lastENode = len(sequences)
allnodes = []

def readTree(s):
	#nonlocal allnodes
	global lastENode
	#nonlocal edges
	if s[0] == '(':
		node = lastENode
		lastENode += 1

		node1, dis1, s = readTree(s[1:])
		if s[0] != ',':
			raise "error"
		node2, dis2, s = readTree(s[1:])
		if s[0] != ')':
			raise Exception("error in reading ): {}".format(s))
		edges.append((node, node1, dis1))
		edges.append((node, node2, dis2))
		allnodes.append(node)
		s = s[1:]
	else:
		r1 = re.match(r"[a-zA-Z]*([0-9]*)", s)
		if not r1:
			raise "error in name of cell"
		node = int(r1.group(1))-1
		allnodes.append(node)
		s = s[len(r1.group()):]
	if len(s) > 0:
		if s[0] != ':':
			raise "error"
		s = s[1:]
		r1 = re.match(r"([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)", s)
		if not r1:
			raise Exception("error in reading number: {}".format(s))
		dis = float(r1.group())
		return node, dis, s[len(r1.group()):]
	return node, 0, s
	

for l in sciteFile:
	m = re.search("best tree = ([^;]*);", l)
	if m:
		readTree(m.group(1))

allnodes = set(allnodes)
print('{}'.format(' '.join([str(x) for x in allnodes])), file=treeFile)
for (v, u, d) in edges:
	# print("{} -> {}".format(v, u))
	# Order of vertices in our tree is reverse of the order of vertices of SCITE
	print("{}->{} {}".format(u, v, d), file=treeFile)

for n in allnodes:
	clones = []
	if n < len(sequences):
		clones.append(n+1)
	print("{} {}".format(n, ' '.join([str(x) for x in clones])), file=cloneFile)

