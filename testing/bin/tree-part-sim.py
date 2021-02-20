import sys, itertools, re, itertools

treeFile = open(sys.argv[1])
cloneFile = open(sys.argv[2])

treeNodes = re.sub(' +', ' ', treeFile.readline().strip()).split()
n = len(treeNodes)

edges = {i:[] for i in treeNodes}
for l in treeFile:
	x = re.sub(' +', ' ', l.strip()).split(' ')
	e, w = x[0], float(x[1])
	(v, u) = e.split('->')
	# print(v,u,w,x)
	## Here we use the distance reported by the algorithm,
	#dist[v][u] = dist[u][v] = w
	edges[v].append((u,w))
	edges[u].append((v,w))
	# dist[v][u] = ???

nodeClone = {}
for l in cloneFile:
	x = re.sub(' +', ' ', l.strip()).split(' ')
	if len(x) == 0: continue
	nodeClone[x[0]] = x[1:]

child = {}
subtree = {}
def dfs(v, p):
	child[v] = []
	subtree[v] = [v]
	for u, w in edges[v]:
		if u != p:
			child[v].append((u,w))
			dfs(u, v)
			subtree[v].extend(subtree[u])

root = treeNodes[0]
dfs(root, -1)

#print(subtree)

allCells = set.union(*[set(nodeClone[v]) for v in treeNodes])
partitions = []
def addPartition(cells):
	cellsComp = sorted(list(allCells - set(cells)))
	cells = sorted(cells)
	#print("* {} : {} // {} {}".format(cells, cellsComp, len(cells), len(cellsComp)))
	if len(cellsComp) < len(cells) or (len(cellsComp) == len(cells) and list(cellsComp) < list(cells)):
		cells, cellsComp = cellsComp, cells
	#print("- {} : {} // {} {}".format(cells, cellsComp, len(cells), len(cellsComp)))
	#print('{}'.format(' '.join([c for c in cells])))
	if len(cells) != 0:
		partitions.append(cells)
	

for v in treeNodes:
	addPartition(list(set.union(*[set(nodeClone[u]) for u in subtree[v]])))
	if len(nodeClone[v]) != 1 or len(child[v]) > 0:
		for cell in nodeClone[v]:
			addPartition([cell])

#print(([p for p in partitions if len(p) > 1]), file=sys.stderr)

partitions.sort()
partitions = list([k for k,_ in itertools.groupby(partitions)])
for part in partitions:
	print('{}'.format(' '.join([c for c in part])))

