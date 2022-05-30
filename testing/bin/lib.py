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
	#dot.graph_attr['spline'] = 'false'
	dot.graph_attr['splines'] = 'compound'
	#dot.graph_attr['splines'] = 'line'
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
		dot.node(treeNode, desc, color=col, fillcolor=fillcol, style="filled", fontcolor=fontcol, gradientangle="0", penwidth="4", shape="circle", margin="0", fontname="Arial")

	for v, u, w in edges:
		tup = treeEdgeLabel(v, u, w)
		if isinstance(tup, tuple):
			label, edgecol = tup[0], tup[1]
		else:
			label, edgecol = tup, "black"
		dot.edge(u, v, label = label, len = str(w), color=edgecol, fontname="Arial")
			#weight=str(w), 

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

def fillTreeNodeMutations(treeNodeCells, treeChildren, mutationInfo, sequences, cellNames, treeRoot, mutIndex = None, type='common:positive:ignore-leaf'):
	allMutationCount = [ [ sum([1 for s in sequences if s[i] == mut]) for mut in [0, 1] ] for i,m in enumerate(mutationInfo) ]
	# print(allMutationCount)

	def dfs(v):
		myCellsStar = []
		if v in treeNodeCells:
			myCellsStar += [int(u)-1 for u in treeNodeCells[v]]
		if v in treeChildren:
			for u,w in treeChildren[v]:
				myCellsStar = dfs(u) + myCellsStar

		nodeMutations[v] = []
		#treeNodeMutations[v] = []
		for i, m in enumerate(mutationInfo):
			subTreeNormal, subTreeMutated = 0, 0
			for c in myCellsStar:
			#if v in treeNodeCells:
			#	for u in treeNodeCells[v]:
				#c = int(u) - 1
				if sequences[c][i] == 0:
					subTreeNormal += 1
				if sequences[c][i] == 1:
					subTreeMutated += 1
			
			#if subTreeMutated > subTreeNormal:
			if subTreeMutated > 0:
				ignoremut = False
				if v in treeChildren:
					for u,w in treeChildren[v]:
						smi = [submut for ii, submut, desc in nodeMutations[u] if i == ii]
						if len(smi) == 1 and smi[0] == subTreeMutated:
							ignoremut = True
				if not ignoremut:
					nodeMutations[v].append((i, subTreeMutated, str(subTreeNormal) +',' + str(len(myCellsStar)) ))
				#treeNodeMutations[v].append(m['gene']+'/'+str(subTreeMutated))
				print("Mut:{} node: {} gene: {} cells:{} subtree: {},{}/{},{} ".format(i, v, m['gene'], ','.join([cellNames[int(c)-1] for c in treeNodeCells[v]]) , subTreeNormal, subTreeMutated, allMutationCount[i][0], allMutationCount[i][1]))
				
			#oddsratio, pvalue = stats.fisher_exact([[subTreeMutated, allMutationCount[i][1] - subTreeMutated], [subTreeNormal, allMutationCount[i][0]-subTreeNormal]], alternative='greater')
			#if pvalue <= 0.01 and v != treeRoot:
			#	treeNodeMutations[v].append(m['gene'])
			#	print([[subTreeMutated, allMutationCount[i][1] - subTreeMutated], [subTreeNormal, allMutationCount[i][0]-subTreeNormal+100]])
			#	print("testing for mutation:{} location: {} p-value: {} gene: {} cells:{} subtree: {},{}/{},{} ".format(i, v, pvalue, m['gene'], ','.join(treeNodeCells[v]) , subTreeNormal, subTreeMutated, allMutationCount[i][0], allMutationCount[i][1]))
		return myCellsStar

	def numberOfMutationsInNodeCells(v, i):
		cellsNormal = cellsMutated = 0
		for c in (treeNodeCells[v] if v in treeNodeCells else []):
			c = int(c)-1
			if sequences[c][i] == 0:
				cellsNormal += 1
			if sequences[c][i] == 1:
				cellsMutated += 1
		return cellsNormal, cellsMutated
	
	def numberOfMutationsInChildren(v, i):
		mutationInChildren = noMutationInChildren = 0
		for u,w in treeChildren[v] if v in treeChildren else []:
			if i in nodeMutationDict[u]: mutationInChildren += 1
			else: noMutationInChildren += 1
		return noMutationInChildren, mutationInChildren
		 

	# Prints all mutations positive and negative on internal nodes. 
	# Caucluation of mutations is somehow greedy.
	def dfsGreedyMinusBias(v):
		if v in treeChildren:
			for u,w in treeChildren[v]:
				dfsGreedyMinusBias(u)

		nodeMutationDict[v] = set()

		nodeMutations[v] = []
		#treeNodeMutations[v] = []
		for i, m in enumerate(mutationInfo):
			cellsNormal, cellsMutated = numberOfMutationsInNodeCells(v, i)
			noMutationInChildren, mutationInChildren = numberOfMutationsInChildren(v, i)
	
			nodeHasMutation = False

			# we have at least one cell and all the cells in this node agree on this mutation
			if (v in treeNodeCells and len(treeNodeCells[v]) > 0 and (cellsNormal == 0 or cellsMutated == 0)):
				if cellsMutated > 0:
					nodeHasMutation = True
			elif (mutationInChildren == 0 or noMutationInChildren == 0) and v in treeChildren and len(treeChildren[v]) > 0:
				#either no cell in this node or the cells do not agree but children agree 
				if mutationInChildren > 0:
					nodeHasMutation = True
			else:
				#niether cells nor all the children agree on a mutation
				if mutationInChildren > noMutationInChildren:
					nodeHasMutation = True
			if nodeHasMutation:
				nodeMutationDict[v].add(i)

			for u,w in treeChildren[v] if v in treeChildren else []:
				if (i in nodeMutationDict[v]) ^ (i in nodeMutationDict[u]):
					nodeMutations[u].append((i, '+' if (i in nodeMutationDict[u]) else '-', ''))

			if v == treeRoot:
				if (i in nodeMutationDict[v]):
					nodeMutations[v].append((i, '+', ''))
	
	def vectorSum(a):
		if len(a) == 0: return ()
		ret = [0 for i in range(len(a[0]))]
		for x in a:
			for i, y in enumerate(x):
				ret[i] += y
		return tuple(ret)


						

	# Prints all mutations common in the whole subtree
	def dfsOnlyCommonMutations(v):
		if v in treeChildren:
			for u,w in treeChildren[v]:
				dfsOnlyCommonMutations(u)

		nodeMutationDict[v] = set()

		nodeMutations[v] = []
		#treeNodeMutations[v] = []
		for i, m in enumerate(mutationInfo):
			cellsNormal, cellsMutated = numberOfMutationsInNodeCells(v, i)
			nodeCellCount = len(treeNodeCells[v])
			# print('d', v, i, cellsNormal, cellsMutated)

			# if cellsMutated > 0:
			# 	nodeMutationDict[v].add((i, 1))
			# if cellsNormal > 0:
			# 	nodeMutationDict[v].add((i, -1))
			nodeCellsCombatibility = {1: cellsNormal, -1: cellsMutated}
			
			# (i,2d) means compatibility of all children of v with d
			for d in [-1, 1]:
				# if (v not in treeChildren or all([(i,2*d) in nodeMutationDict[u] for u,w in treeChildren[v]])) and \
				# 	(i,-d) not in nodeMutationDict[v]:
				# 	nodeMutationDict[v].add((i, 2*d))
				if v in treeChildren:
					r = vectorSum([(nmd[2], nmd[3]) for u,w in treeChildren[v] for nmd in nodeMutationDict[u] if nmd[0] == i and nmd[1] == 2*d])
				else:
					r = (0,0)
				# compatible and non missing vaules for this node
				r = vectorSum([r, (nodeCellCount-nodeCellsCombatibility[-d], nodeCellCount)])
				nodeMutationDict[v].add((i, 2*d, r[0], r[1]))

			def showMutationOnChild(v, u, d):
				if v != -1:
					nmdv = [(nmd[2], nmd[3]) for nmd in nodeMutationDict[v] if nmd[0] == i and nmd[1] == 2*d][0]
				else:
					nmdv = (0, 1)
				nmdu = [(nmd[2], nmd[3]) for nmd in nodeMutationDict[u] if nmd[0] == i and nmd[1] == 2*d][0]
				# Note: if nmdv[1] == 0, for its children we cannot have nmdu[1] > 0!
				return nmdv[0] < threshold * nmdv[1] and nmdu[0] >= threshold * nmdu[1]

			for d in [-1, 1]:
				for u,w in treeChildren[v] if v in treeChildren else []:
					# if ((i,2*d) in nodeMutationDict[v]) ^ ((i,2*d) in nodeMutationDict[u]) and\
					# 	((i,2*d) not in nodeMutationDict[u] or (i,-2*d) not in nodeMutationDict[u]) and\
					if showMutationOnChild(v,u,d) and not showMutationOnChild(v,u,-d) and\
						(('ignore-leaf' not in commandRest) or (u in treeChildren and len(treeChildren[u])>0)):
						if commandSubtype == 'all' or (d == 1 and commandSubtype == 'positive'):
							nodeMutations[u].append((i, '+' if d == 1 else '-', ''))
						if commandSubtype == 'all' or (d == -1 and commandSubtype == 'negative'):
							nodeMutations[u].append((i, '+' if d == 1 else '-', ''))

			if v == treeRoot:
				for d in [-1, 1]:
					# if ((i,2*d) in nodeMutationDict[v]):
					if showMutationOnChild(-1,v,d) and not showMutationOnChild(-1,v,-d):
						if commandSubtype == 'all' or (d == 1 and commandSubtype == 'positive'):
							nodeMutations[v].append((i, '+' if d == 1 else '-', ''))
						if commandSubtype == 'all' or (d == -1 and commandSubtype == 'negative'):
							nodeMutations[v].append((i, '+' if d == 1 else '-', ''))
		# print('d', v, nodeMutationDict[v])


	# Prints all mutations common in the node and its direct children. If a node does not have a sample, we include recursively its children as direct children of the parent
	def dfsOnlyCommonMutationsInDirectChildren(v):
		def directChildrenWithSample(v, canStopHere = False):
			ret = []
			if canStopHere:
				ret.append(v)
				if len(treeNodeCells[v]) > 0: return ret
			if v in treeChildren:
				for u,w in treeChildren[v]:
					ret += directChildrenWithSample(u, canStopHere = True)
			return ret

		if v in treeChildren:
			for u,w in treeChildren[v]:
				dfsOnlyCommonMutationsInDirectChildren(u)

		nodeMutationDict[v] = set()

		nodeMutations[v] = []
		#treeNodeMutations[v] = []
		for i, m in enumerate(mutationInfo):
			cellsNormal, cellsMutated = numberOfMutationsInNodeCells(v, i)
			nodeCellCount = len(treeNodeCells[v]) if v in treeNodeCells else 0
			# print('d', v, i, cellsNormal, cellsMutated)

			# if cellsMutated > 0:
			# 	nodeMutationDict[v].add((i, 1))
			# if cellsNormal > 0:
			# 	nodeMutationDict[v].add((i, -1))
			if cellsNormal == 0:
				nodeMutationDict[v].add((i, -1))
			if cellsMutated == 0:
				nodeMutationDict[v].add((i, 1))
			nodeCellsCombatibility = {1: cellsNormal, -1: cellsMutated}
			
			# (i,2d) means compatibility of all children of v with d
			for d in [-1, 1]:
				# if (v not in treeChildren or all([(i,2*d) in nodeMutationDict[u] for u,w in treeChildren[v]])) and \
				# 	(i,-d) not in nodeMutationDict[v]:
				# 	nodeMutationDict[v].add((i, 2*d))
				if v in treeChildren:
					r = all([ (i,d) in nodeMutationDict[u] for u in directChildrenWithSample(v)])
				else:
					r = True
				# compatible and non missing vaules for this node
				if r:
					nodeMutationDict[v].add((i, 2*d))
			
			for d in [-1, 1]:
				if ((i,2*d) in nodeMutationDict[v] and (i,-2*d) not in nodeMutationDict[v]) and \
					(('ignore-leaf' not in commandRest) or (v in treeChildren and len(treeChildren[v])>0)):
					if commandSubtype == 'all' or (d == 1 and commandSubtype == 'positive'):
						nodeMutations[v].append((i, '+' if d == 1 else '-', ''))
					if commandSubtype == 'all' or (d == -1 and commandSubtype == 'negative'):
						nodeMutations[v].append((i, '+' if d == 1 else '-', ''))

	nodeMutations = {}
	nodeMutationDict = {}
	#dfs(treeRoot)
	# dfsGreedyMinusBias(treeRoot)
	commandType, commandSubtype, commandRest = type.split(':')
	if commandType.startswith('common='):
		threshold = float(commandType.split('=')[1])
		dfsOnlyCommonMutations(treeRoot)
	elif commandType == 'direct-children':
		# threshold = float(commandType.split('=')[1])
		dfsOnlyCommonMutationsInDirectChildren(treeRoot)
	else:
		raise RuntimeError("Invalid command type: " + commandType)
	treeNodeMutations = {}
	for v, muts in nodeMutations.items():
		treeNodeMutations[v] = [ (mutationInfo[i]['gene'] if mutationInfo[i]['gene'] != '-' else 'M#'+str(i+1)) + '' + str(mut) + (',' + str(desc) if len(desc) > 0 else '') for i, mut, desc in muts if mutIndex is None or i == mutIndex]
	return treeNodeMutations
