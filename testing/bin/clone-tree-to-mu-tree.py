import sys, re, argparse
from graphviz import Digraph
from lib import *
import scipy.stats as stats

parser = argparse.ArgumentParser(description='Clone tree to mutation tree.')
parser.add_argument('treeFileName')
parser.add_argument('cloneFileName')
parser.add_argument('seqFileName')
parser.add_argument('mutationInfoFileName')
parser.add_argument('cellNamesFileName')
parser.add_argument('outputFileName')
parser.add_argument('--compress', dest='compress', action='store_const', const = True, default = False)
parser.add_argument('--mark-mutations', dest='markMutations', action='store_const', const = True, default = False)
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
if args.markMutations:
        allMutationCount = [ [ sum([1 for s in sequences if s[i] == mut]) for mut in [0, 1] ] for i,m in enumerate(mutationInfo) ]
        # print(allMutationCount)

        def dfs(v):
                myCellsStar = []
                if v in treeNodeCells:
                        myCellsStar += [int(u)-1 for u in treeNodeCells[v]]
                if v in treeChildren:
                        for u,w in treeChildren[v]:
                                myCellsStar = dfs(u) + myCellsStar

                treeNodeMutations[v] = []
                for i, m in enumerate(mutationInfo):
                        subTreeNormal, subTreeMutated = 0, 0
                        for c in myCellsStar:
                                if sequences[c][i] == 0:
                                        subTreeNormal += 1
                                if sequences[c][i] == 1:
                                        subTreeMutated += 1
                        
                        oddsratio, pvalue = stats.fisher_exact([[subTreeMutated, allMutationCount[i][1] - subTreeMutated], [subTreeNormal, allMutationCount[i][0]-subTreeNormal]], alternative='greater')
                        # if i == 279:
                        #         print([[subTreeMutated, allMutationCount[i][1] - subTreeMutated], [subTreeNormal, allMutationCount[i][0]-subTreeNormal]])
                        #         print("testing for mutation:{} location: {} p-value: {} gene: {} cells:{} subtree: {},{}/{},{} ".format(i, v, pvalue, m['gene'], ','.join(treeNodeCells[v]) , subTreeNormal, subTreeMutated, allMutationCount[i][0], allMutationCount[i][1]))
                        if pvalue <= 0.01 and v != treeRoot:
                                treeNodeMutations[v].append(m['gene'])
                                print([[subTreeMutated, allMutationCount[i][1] - subTreeMutated], [subTreeNormal, allMutationCount[i][0]-subTreeNormal+100]])
                                print("testing for mutation:{} location: {} p-value: {} gene: {} cells:{} subtree: {},{}/{},{} ".format(i, v, pvalue, m['gene'], ','.join(treeNodeCells[v]) , subTreeNormal, subTreeMutated, allMutationCount[i][0], allMutationCount[i][1]))
                return myCellsStar

        dfs(treeRoot)


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
        desc = ', '.join([cellNames[int(c)-1] for c in cells]) 
        if treeNodeMutations is not None and treeNode in treeNodeMutations:
              mutations = sorted(list(set(treeNodeMutations[treeNode])))
              if len(mutations) > 0:
                desc += " ["
                desc += " ".join(mutations)
                desc += "]"
        col="black"
        return desc, col

def treeEdgeLabel(v, u, w):
        return "{:.1f}".format(float(w))

writeGraph(args.outputFileName, treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel)

