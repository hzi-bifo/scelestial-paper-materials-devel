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
parser.add_argument('--mark-mutations', dest='markMutations', help='example: common=1:positive:ignore-leaf')
parser.add_argument('--mark-mutations-separated', dest='markMutationsSeparated', action='store_const', const = True, default = False)
parser.add_argument('--merge-colons', dest='margeClones', action='store_const', const = True, default = False)
parser.add_argument('--cell-type', dest='cellTypeFileName')

parser.add_argument('--compressed', )
args = parser.parse_args()

vertices, edges, treeParent, treeChildren, treeRoot = loadTree(args.treeFileName)
treeNodeCells = loadClones(args.cloneFileName)

sequences = loadSequenceFile(args.seqFileName)
mutationInfo = loadMutationInfoFile(args.mutationInfoFileName)
cellNames = loadCellNames(args.cellNamesFileName)

cellTypes = None
if args.cellTypeFileName:
	cellTypes = loadCellTypes(args.cellTypeFileName)

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
if args.markMutations:
	# treeNodeMutations = fillTreeNodeMutations()
	treeNodeMutations = fillTreeNodeMutations(treeNodeCells, treeChildren, mutationInfo, sequences, cellNames, treeRoot, None, args.markMutations)

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
		#print("{}/{}".format(cells, len(cellNames)), file=sys.stderr)
		if cellTypes is not None:
			CN = set([cellTypes[int(c)-1] == 0 for c in cells])
		else:
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
			desc += x[0] + x[1]
			r += len(x[0] + x[1])

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
		seqmut = math.sqrt(len(mutations)/6)
		if len(mutations) > 0:
			desc += " ["
			lastLineLen = 0
			for i, mut in enumerate(mutations):
				if lastLineLen >= seqmut:
					desc += '\n'
					lastLineLen = 0
				desc += mut + " "
				lastLineLen += 1
			#desc += " ".join(mutations)
			desc += "]"
	return desc, col, fontcol, fillcol

def treeEdgeLabel(v, u, w):
	return "", "#A2A2A2"
#return "{:.1f}".format(float(w))

if args.markMutationsSeparated:
	for i, mut in enumerate(mutationInfo):
		# treeNodeMutations = fillTreeNodeMutations(i)
		treeNodeMutations = fillTreeNodeMutations(treeNodeCells, treeChildren, mutationInfo, sequences, cellNames, treeRoot, i, args.markMutations)
		writeGraph(args.outputFileName + '-' + str(i), treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel)
else:
	writeGraph(args.outputFileName, treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel)

