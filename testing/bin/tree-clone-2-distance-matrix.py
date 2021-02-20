import sys, itertools, re

treeFile = open(sys.argv[1])
cloneFile = open(sys.argv[2])

normalize = 'no'
argc = 3
while argc < len(sys.argv):
	if sys.argv[argc] == '-no-normal':
		normalize = 'no'
	elif sys.argv[argc] == '-normal-mean':
		normalize = 'mean'
	argc += 1

names = re.sub(' +', ' ', treeFile.readline().strip()).split()
n = len(names)

INF = 1e10
dist = { i:{ j:INF for j in names} for i in names }
for i in names:
	# for j in names:
	#         dist[i][i] = INF
	dist[i][i] = 0

edges = {i:[] for i in names}
dedges = {i:[] for i in names}
for l in treeFile:
	x = re.sub(' +', ' ', l.strip()).split(' ')
	e, w = x[0], float(x[1])
	(v, u) = e.split('->')
	# print(v,u,w,x)
	## Here we use the distance reported by the algorithm,
	dist[v][u] = dist[u][v] = w
	edges[v].append((u,w))
	edges[u].append((v,w))
	dedges[v].append((u,w))
	## I changed it to the following:
	# dist[v][u] = ???

#print('Floyd ... n={}'.format(len(names)), file=sys.stderr)
#for k in names:
#        for i in names:
#                for j in names:
#                        dist[i][j] = min( dist[i][j] , dist[i][k] + dist[k][j] )
#print('Floyd done', file=sys.stderr)

clone = {}
for l in cloneFile:
	x = re.sub(' +', ' ', l.strip()).split(' ')
	if len(x) == 0: continue
	for se in x[1:]:
		clone[se] = x[0]

# print(clone)


inDeg = {i:0 for i in names}
for u, eds in dedges.items():
	for ed in eds:
		v, w = ed
		inDeg[v] += 1
roots = [v for v in names if inDeg[v] == 0]

def dfs(i, mark, v, d):
	mark[v] = True
	dist[i][v] = d
	for u, w in edges[v]:
		if u not in mark:
			dfs(i, mark, u, d+w)

if normalize == 'mean':
	if len(roots) != 1:
		print('More than one root!: {}'.format(roots), file=sys.stderr)
	root = roots[0]
	mark = {}
	dfs(root, mark, root, 0)
	rootToCellsDistance = sum([dist[root][c] for cell, c in clone.items()])
	newEdges = {}
	for u, eds in edges.items():
		neds = [(ed[0], ed[1]/rootToCellsDistance) for ed in eds]
		newEdges[u] = neds
	edges = newEdges


for i in names:
	mark = {}
	dfs(i, mark, i, 0)

# print(dist)

columns = [str(k) for k in sorted([int(k) for k in clone.keys()]) ]
# print(columns, file=sys.stderr)


print('""', end=" ")
for col in columns:
	print('"{}"'.format(col), end=" ")
print()
for col in columns:
	print('"{}"'.format(col), end=" ")
	for j in columns:
		print(dist[clone[col]][clone[j]], end=" ")
	print()

# print('""', end=" ")
# for (j,d) in clone.items():
#         print('"{}"'.format(j), end=" ")
# print()
# for (i,c) in clone.items():
#         print('"{}"'.format(i), end=" ")
#         for (j,d) in clone.items():
#                 print(dist[c][d], end=" ")
#         print()
