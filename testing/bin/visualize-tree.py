import sys, re, argparse
from graphviz import Digraph
from lib import *
import math
import importlib  


parser = argparse.ArgumentParser(description='Generate samples from single-cell data set.')
parser.add_argument('title')
parser.add_argument('treeFileName')
parser.add_argument('cloneFileName')
parser.add_argument('outputFileName')
parser.add_argument('--type-file', dest='typeFile', action='store', required = False)
parser.add_argument('--compressed-out-file', dest='compressedOutFileName', action='store', required = False)
parser.add_argument('--cell-name', dest='cellNamesFileName')
parser.add_argument('--color', dest='colorFileName')
parser.add_argument('--seq', dest='seqFileName', action='store', required=False)
parser.add_argument('--mark-mutations', dest='markMutations', help='example: common=1:positive:ignore-leaf')
parser.add_argument('--mutation-info', dest='mutationInfoFileName', action='store', required=False)
args = parser.parse_args()

title = args.title
classColor = ["aquamarine", "azure", "beige", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "cyan", "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgreen", "darkgrey", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey", "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite", "forestgreen", "fuchsia", "gainsboro", "ghostwhite", "gold", "goldenrod", "gray", "grey", "green", "greenyellow", "honeydew", "hotpink", "indianred", "indigo", "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightgray", "lightgreen", "lightgrey", "lightpink", "lightsalmon", "lightseagreen", "lightskyblue", "lightslategray", "lightslategrey", "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "magenta", "maroon", "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab", "orange", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple", "red", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "silver", "skyblue", "slateblue", "slategray", "slategrey", "snow", "springgreen", "steelblue", "tan", "teal", "thistle", "tomato", "turquoise", "violet", "wheat", "white", "whitesmoke", "yellow", "yellowgreen"]
classColor = {i:s for i, s in enumerate(classColor)}
classColor['general-col'], classColor['general-fill'] = 'gray', 'white'
classColor['mix-col'], classColor['mix-fill'] = 'gray', 'white'

vertices, edges, treeParent, treeChildren, treeRoot = loadTree(args.treeFileName)
treeNodeCells = loadClones(args.cloneFileName)
cellNames = loadCellNames(args.cellNamesFileName) if args.cellNamesFileName else None
if args.typeFile != None:
	headerClasses = loadCellTypes(args.typeFile)
if args.colorFileName != None:
	classColor_ = loadTable(args.colorFileName)
	#print(classColor_)
	classColor = {int(r[0]) if r[0].isnumeric() else r[0]:r[1].strip() for r in classColor_}
if args.seqFileName:
	sequences = loadSequenceFile(args.seqFileName)
treeNodeMutations = {}
if args.mutationInfoFileName:
	mutationInfo = loadMutationInfoFile(args.mutationInfoFileName)
if args.markMutations:
	if 'mutationInfo' not in vars():
		mutationInfo = []
		for i in range(len(sequences[0])):
			mutationInfo.append({'gene' : '-', 'geneInfo' : '', id: str(i+1)})
	treeNodeMutations = fillTreeNodeMutations(treeNodeCells, treeChildren, mutationInfo, sequences, cellNames, treeRoot, None, args.markMutations)


def treeNodeDescColor(treeNode, cells):
	desc = ','.join(cells)
	#col, fillcol = "grey", "white"
	col, fillcol = classColor['general-col'], classColor['general-fill']
	if args.typeFile != None:
		for cell in cells[0:]:
			if int(cell)-1 >= len(headerClasses):
				print('Invalid cell classe, cell={}'.format(cell), file=sys.stderr)
				print(headerClasses, file=sys.stderr)
		desc = ','.join(["{}/{}:{}".format(treeNode, cell, headerClasses[int(cell)-1]) for cell in cells[0:]])
		#if len(cells) == 1:
		if len(cells) == 0:
			#col, fillcol = "grey", "white"
			col, fillcol = classColor['general-col'], classColor['general-fill']
		if len(set([headerClasses[int(c)-1] for c in cells])) == 1:
			fillcol = col = classColor[headerClasses[int(cells[0])-1]]
		if len(set([headerClasses[int(c)-1] for c in cells])) > 1:
			col, fillcol = classColor['mix-col'], classColor['mix-fill']
		w = 6 * len(cells)
		desc = ''
		r = 0
		for c in cells:
			x = cellNames[int(c)-1] if cellNames is not None else c
			if len(desc) > 0: 
				if r >= math.sqrt(w):
					desc += ',\n'
					r = 0
				else:
					desc += ', '
					r += 2
			desc += x
			r += len(x)
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
	return desc, col, "black", fillcol

def treeEdgeLabel(v, u, w):
	return "{:.1f}".format(float(w))

writeGraph(args.outputFileName, treeNodeCells, list(treeNodeCells.keys()), edges, treeNodeDescColor, treeEdgeLabel)

if args.compressedOutFileName is not None:
	dot = Digraph(comment=title, format='pdf')
	dot.graph_attr['rankdir'] = 'LR'

	print('{}'.format(treeRoot), file=sys.stderr)

	compressedTreeChildren = compressedTree(treeChildren, treeNodeCells, treeRoot)
	for v, children in compressedTreeChildren.items():
		desc = ""
		if v in treeNodeCells:
			desc = ','.join(treeNodeCells[v])
		dot.node(v, desc, color="black")
		for c, cw in children:
			# print(" adding edge {} {} {}".format(type(c), type(v), type(cw)))
			dot.edge(c, v, weight=str(cw), label = "{:.1f}".format(cw))

	dot.render(args.compressedOutFileName)

