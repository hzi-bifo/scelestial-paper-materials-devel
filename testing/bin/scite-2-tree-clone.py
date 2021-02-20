import sys, itertools, re

# def distance(sequences, s1, s2):
#         if s1 >= len(sequences):
#                 a = [0 for i in range(len(sequences[0]))]
#         else:
#                 a = sequences[s1]
#         if s2 >= len(sequences):
#                 b = [0 for i in range(len(sequences[0]))]
#         else:
#                 b = sequences[s2]
#         d = 0
#         for x, y in zip(a, b):
#                 if x != y:
#                         d += 1
#         return d

def distance(seq, mutationSet):
        # print(len(seq))
        # print(seq)
        # print(mutationSet)
        s2 = [0 for s in seq]
        for m in mutationSet:
                if m-1 < len(seq):
                        s2[m-1] = 1
        d = 0
        for x, y in zip(seq, s2):
                if x != y and x != 3:
                        d += 1
        return d

sciteFile = open(sys.argv[1])
seqFile = open(sys.argv[2])
treeFile = open(sys.argv[3], "w+")
cloneFile = open(sys.argv[4], "w+")

sequences = []
for l in seqFile:
        for num, val in enumerate(l.strip().split()):
                while num >= len(sequences):
                        sequences.append([])
                sequences[num].append(int(val))

muts = []
edges = []
parent = {}
for l in sciteFile:
        m = re.search("(\\d+) -> (\\d+);", l)
        if m:
                v, u = int(m.group(1)), int(m.group(2))
                edges.append((v,u))
                muts.append(v)
                muts.append(u)
                parent[u] = v

muts = set(muts)
for c in muts:
        print(c, file=treeFile, end = " ")
print(file = treeFile)

for (v, u) in edges:
        # print("{} -> {}".format(v, u))
        # Order of vertices in our tree is reverse of the order of vertices of SCITE
        print("{}->{} {}".format(u, v, 1), file=treeFile)

allMutations = {}
mark = {}

def dfs(m):
        mark[m] = True
        if m not in parent:
                allMutations[m] = [m]
                return
        if parent[m] not in mark:
                dfs(parent[m])
        allMutations[m] = allMutations[parent[m]] + [m]

for m in muts:
        if m not in mark:
                dfs(m)

clones = { m:[] for m in muts }

for s, seq in enumerate(sequences):
        minD = len(muts) * 2 + 10
        for m in muts:
                d = distance(seq, allMutations[m])
                if d < minD:
                        minD = d
                        minM = m
        clones[minM].append(s)

for cl, cells in clones.items():
        print("{0}".format(cl), end=" ", file=cloneFile)
        for c in cells:
                print("{0}".format(c+1), end = " ", file=cloneFile)
        print(file=cloneFile)

