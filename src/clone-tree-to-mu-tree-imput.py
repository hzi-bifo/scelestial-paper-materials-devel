import sys, re, argparse
from graphviz import Digraph
from lib import *
import scipy.stats as stats
import math

BLUE = "#69c5f0"
BLACK = "black"
ORANGE = "#f57433"
BROWN = "#af9d92"

parser = argparse.ArgumentParser(description='Clone tree to mutation tree.')
parser.add_argument('treeFileName')
parser.add_argument('cloneFileName')
parser.add_argument('seqFileName')
parser.add_argument('mutationInfoFileName')
parser.add_argument('cellNamesFileName')
parser.add_argument('outputFileName')
parser.add_argument('--compress', dest='compress', action='store_const', const = True, default = False)
parser.add_argument('--mark-mutations', dest='markMutations', action='store_const', const = True, default = False)
parser.add_argument('--mark-mutations-separated', dest='markMutationsSeparated', action='store_const', const = True, default = False)
parser.add_argument('--merge-colons', dest='margeClones', action='store_const', const = True, default = False)

parser.add_argument('--compressed', )
args = parser.parse_args()

vertices, edges, treeParent, treeChildren, treeRoot = loadTree(args.treeFileName)
treeNodeCells = loadClones(args.cloneFileName)

sequences = loadSequenceFile(args.seqFileName)
mutationInfo = loadMutationInfoFile(args.mutationInfoFileName)
cellNames = loadCellNames(args.cellNamesFileName)

if args.margeClones:
	## May be better to be after compressing the tree

	def sameColon(v, u, w):
		m = len(sequences[0])
		pv = stats.binom_test(w, m, 0.2, alternative='less')
		print("same colon test: {}/{} = {}".format(w, m, pv))
		return pv < 0.01
	mergeParent = { v:v for v in treeChildren.keys() }
	toBeMerged = {}
	for v, childrenDistList in treeChildren.items():
		for c, w in childrenDistList:
			if sameColon(v, c, w):
				mergeParent[c] = v
				toBeMerged[c] = True
	
	newTreeChildren = { v : [] for v in treeChildren.keys() }
	newTreeNodeCells = { v : [] for v in treeNodeCells.keys() }
	def dfs(v, firstKeptNode, wToHere):
		newTreeNodeCells[firstKeptNode] += treeNodeCells[v]
		if v in treeChildren:
			for u, w in treeChildren[v]:
				if u not in toBeMerged:
					newTreeChildren[firstKeptNode].append((u, wToHere + w))
					dfs(u, u, 0)
				else:
					dfs(u, firstKeptNode, wToHere + w)
	dfs(treeRoot, treeRoot, 0)

	for v, childrenDistList in treeChildren.items():
		while mergeParent[mergeParent[v]] != mergeParent[v]:
			mergeParent[v] = mergeParent[mergeParent[v]]
		newTreeChildren[v] 


	treeNodeCells = newTreeNodeCells
	treeChildren = newTreeChildren
	treeRoot = treeRoot

treeNodeMutations = {}
def fillTreeNodeMutations(mutIndex = None):
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

	nodeMutations = {}
	dfs(treeRoot)
	for v, muts in nodeMutations.items():
		treeNodeMutations[v] = [ mutationInfo[i]['gene'] + '/' + str(mut) + ',' + str(desc) for i, mut, desc in muts if mutIndex is None or i == mutIndex]

if args.markMutations:
	fillTreeNodeMutations()

nodes = list(treeNodeCells.keys())
if args.compress:
	compressedTreeChildren = compressedTree(treeChildren, treeNodeCells, treeRoot)
	edges = []
	for par, cwList in compressedTreeChildren.items():
		for c, w in cwList:
			# print("{} -> {}".format(par, cw))
			edges.append((c, par, w))
	nodes = list(compressedTreeChildren.keys())


def treeNodeDescColor(treeNode, cells):
	#desc = ', '.join([cellNames[int(c)-1] for c in cells]) 
	desc = ''
	col, fontcol, fillcol = "#C2C2C2", "black", "white"
	if len(cells) > 0:
		CN = set([cellNames[int(c)-1].find('BC') != -1 for c in cells])
		#desc = ', '.join([cellNames[int(c)-1].split('-')[0][1] + cellNames[int(c)-1].split('-')[1] for c in cells])
		w = 4 * len(cells)
		desc = ''
		r = 0
		for c in cells:
			x = cellNames[int(c)-1].split('-')
			if len(desc) > 0: 
				if r >= math.sqrt(w):
					desc += ',\n'
					r = 0
				else:
					desc += ', '
					r += 2
			desc += x[0][1] + x[1]
			r += len(x[0][1] + x[1])

		if len(CN) == 1:
			cT = CN.pop()
			if cT == False:
				col, fontcol, fillcol = BLUE, BLACK, BLUE
			else:
				col, fontcol, fillcol = ORANGE, "black", ORANGE
		else:
			col, fontcol, fillcol = BROWN, "black", BROWN
	if treeNodeMutations is not None and treeNode in treeNodeMutations:
		mutations = sorted(list(set(treeNodeMutations[treeNode])))
		seqmut = math.sqrt(len(mutations)/7)
		if len(mutations) > 0:
			desc += " ["
			lastLineLen = 0
			for i, mut in enumerate(mutations):
				desc += mut + " "
				lastLineLen += 1
				if lastLineLen >= seqmut:
					desc += '\n'
					lastLineLen = 0
			#desc += " ".join(mutations)
			desc += "]"
	return desc, col, fontcol, fillcol

def treeEdgeLabel(v, u, w):
	return "", "#A2A2A2"
#return "{:.1f}".format(float(w))

if args.markMutationsSeparated:
	for i, mut in enumerate(mutationInfo):
		fillTreeNodeMutations(i)
		writeGraph(args.outputFileName + '-' + str(i), treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel)
else:
	writeGraph(args.outputFileName, treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel)

