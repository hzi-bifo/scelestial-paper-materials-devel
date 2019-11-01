import sys, itertools, re

treeFile = open(sys.argv[1])
cloneFile = open(sys.argv[2])

names = re.sub(' +', ' ', treeFile.readline().strip()).split()
n = len(names)

INF = 1e10
dist = { i:{ j:INF for j in names} for i in names }
for i in names:
	# for j in names:
	#         dist[i][i] = INF
	dist[i][i] = 0

edges = {i:[] for i in names}
for l in treeFile:
	x = re.sub(' +', ' ', l.strip()).split(' ')
	e, w = x[0], float(x[1])
	(v, u) = e.split('->')
	# print(v,u,w,x)
	## Here we use the distance reported by the algorithm,
	dist[v][u] = dist[u][v] = w
	edges[v].append((u,w))
	edges[u].append((v,w))
	## I changed it to the following:
	# dist[v][u] = ???

#print('Floyd ... n={}'.format(len(names)), file=sys.stderr)
#for k in names:
#        for i in names:
#                for j in names:
#                        dist[i][j] = min( dist[i][j] , dist[i][k] + dist[k][j] )
#print('Floyd done', file=sys.stderr)

for i in names:
	mark = {}
	def dfs(v, d):
		mark[v] = True
		dist[i][v] = d
		for u, w in edges[v]:
			if u not in mark:
				dfs(u, d+w)
	dfs(i, 0)

# print(dist)

clone = {}
for l in cloneFile:
	x = re.sub(' +', ' ', l.strip()).split(' ')
	if len(x) == 0: continue
	for se in x[1:]:
		clone[se] = x[0]

# print(clone)

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
