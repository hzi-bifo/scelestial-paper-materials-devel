from graphviz import Digraph

def loadTree(treeFileName):
	treeFile = open(treeFileName)
	l = treeFile.readline()
	vertices = l.strip().split()
	edges = []
	treeParent = {}
	treeChildren = {}
	for l in treeFile:
		(e, w) = l.strip().split(' ')
		(v, u) = e.split('->')
		edges.append((v,u,float(w)))

		treeParent[v] = u
		if u not in treeChildren: treeChildren[u] = []
		treeChildren[u].append((v, float(w)))
	

	treeRoot = list(treeParent.keys())[0]
	while treeRoot in treeParent:
		treeRoot = treeParent[treeRoot]


	return vertices, edges, treeParent, treeChildren, treeRoot

def loadClones(cloneFileName):
	cloneFile = open(cloneFileName)
	treeNodeCells = {}

	for l in cloneFile:
		x = l.strip().split()
		treeNodeCells[x[0]] = x[1:]
	
	return treeNodeCells

def loadCellTypes(cellTypeFileName):
	headerClasses = []
	f = open(cellTypeFileName)
	for l in f:
		headerClasses.append(int(l.strip()))
	return headerClasses

def compressedTree(treeChildren, treeNodeCells, treeRoot):
	compressedTreeChildren = {}

	def dfs(v):
		ret = (None, None)
		children = []
		if v in treeChildren:
			for u, cw in treeChildren[v]:
				c, w = dfs(u)
				if c is not None:
					w += cw
					children.append((c,w))
					ret = (c,w)
		
		if (v in treeNodeCells and len(treeNodeCells[v]) > 0) or len(children) > 1 or v == treeRoot:
			# desc = ""
			# if v in treeNodeCells:
			#	 desc = ','.join(treeNodeCells[v])
			# dot.node(v, desc, color=col)
			ret = (v, 0)

		if ret == (v, 0):
			compressedTreeChildren[v] = []
			for c, cw in children:
				# print(" adding edge {} {} {}".format(type(c), type(v), type(cw)))
				# dot.edge(c, v, weight=str(cw), label = "{:.1f}".format(cw))
			
				compressedTreeChildren[v].append((c, cw))

		return ret
		

	dfs(treeRoot)
	return compressedTreeChildren

## File should be in SCITE format, i.e rows are sites and columns are cells (samples)
## returns a list of cells (samples), sequences[i][j] corresponds to cell i (zero-based), mutation j (zero-based)
def loadSequenceFile(seqFileName):
	sequences = []
	seqFile = open(seqFileName)
	for l in seqFile:
		for num, val in enumerate(l.strip().split()):
			while num >= len(sequences):
				sequences.append([])
			sequences[num].append(int(val))
	return sequences
	
def loadMutationInfoFile(mutationInfoFileName):
	mutationInfoFile = open(mutationInfoFileName)
	mutationInfo = []
	for l in mutationInfoFile:
		x = l.strip().split("\t")
		mutationInfo.append({'gene' : x[11], 'geneInfo' : x[12], id: x[0]})
	return mutationInfo

def writeGraph(outputFileName, treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel):
	dot = Digraph(format='pdf')
	dot.graph_attr['rankdir'] = 'LR'
	for treeNode in nodes:
		# for treeNode, cells in treeNodeCells.items():
		cells = []
		if treeNode in treeNodeCells:
			cells = treeNodeCells[treeNode]
		prop = treeNodeDescColor(treeNode, cells)
		if len(prop) == 2:
			desc, col, fontcol, fillcol = prop[0], prop[1], "black", "none"
		else:
			desc, col, fontcol, fillcol = prop[0], prop[1], prop[2], prop[3]
		dot.node(treeNode, desc, color=col, fillcolor=fillcol, style="filled", fontcolor=fontcol, gradientangle="0", penwidth="4", shape="circle", margin="0")

	for v, u, w in edges:
		tup = treeEdgeLabel(v, u, w)
		if isinstance(tup, tuple):
			label, edgecol = tup[0], tup[1]
		else:
			label, edgecol = tup, "black"
		dot.edge(u, v, weight=str(w), label = label, color=edgecol)

	dot.render(outputFileName)

def loadFileRows(fileName):
	f = open(fileName)
	r = []
	for l in f:
		r.append(l.strip())
	return r

def loadTable(fileName, sep="\t"):
	f = open(fileName)
	r = []
	for l in f:
		r.append(l.split(sep))
	return r

def loadCellNames(cellNamesFileName):
	return loadFileRows(cellNamesFileName)

def writeSequenceFile(sequences, fileName):
	f = open(fileName, "w")
	if len(sequences) > 0:
		for i in range(len(sequences[0])):
			print(" ".join([str(seq[i]) for seq in sequences]), file=f)
	f.close()

