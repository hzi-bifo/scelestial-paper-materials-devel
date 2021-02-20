import sys, re, argparse
from graphviz import Digraph
from lib import *
import scipy.stats as stats
import math

parser = argparse.ArgumentParser(description='my tree to fig tree')
parser.add_argument('treeFileName')
parser.add_argument('cloneFileName')
parser.add_argument('seqFileName')
parser.add_argument('outputFileName')

args = parser.parse_args()

vertices, edges, treeParent, treeChildren, treeRoot = loadTree(args.treeFileName)
treeNodeCells = loadClones(args.cloneFileName)

sequences = loadSequenceFile(args.seqFileName)

def writeFigTree(fName, treeRoot, treeChildren, treeNodeCells):
	f = open(fName, 'w+')
	print("#DATASET", file=f)
	print("Begin trees;", file=f)
	print("\tTranslate", file=f)
	allCells = []
	for v, cells in treeNodeCells.items():
		allCells.extend(cells)
	allCells = set(allCells)
	lines = []
	for c in allCells:
		lines.append("\t\t{} {}".format(c, 'Cell'+c))
	print(',\n'.join(lines), file=f)
	print(";", file=f)


	def printMainTree(v, pw):
		#print("printing ... {} {}".format(v, pw), file=sys.stderr)
		if v not in treeChildren or len(treeChildren[v]) == 0:
			return "{}:{}".format('_'.join(treeNodeCells[str(v)]), pw)
		else:
			r = []
			for (cv, cw) in treeChildren[v]:
				r.append(printMainTree(cv, cw))
			return "(" + ','.join(r) + "):{}".format(pw)
			
	print("tree TREE1 = {};".format(printMainTree(treeRoot, 0.0)), file = f)
	print("End;", file = f)
	f.close()
	

print(treeChildren.keys())
writeFigTree(args.outputFileName, treeRoot, treeChildren, treeNodeCells)
