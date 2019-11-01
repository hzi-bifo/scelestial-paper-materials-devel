import sys, itertools, re

seqSciteFile = open(sys.argv[1])
steinerFile = open(sys.argv[2])
treeFile = open(sys.argv[3], "w+")
cloneFile = open(sys.argv[4], "w+")


##Read Seqeunce file and find the one with max 0z
# sequences = []
# for l in seqSciteFile:
#         for num, val in enumerate(l.strip().split()):
#                 while num >= len(sequences):
#                         sequences.append([])
#                 sequences[num].append(int(val))

# maxCnt = -1
# for (idx, seq) in enumerate(sequences):
#         cnt = 0
#         for x in seq:
#                 if x == 0:
#                         cnt+=1
#         if cnt > maxCnt:
#                 maxCnt = cnt
#                 maxSeq = seq
#                 maxSeqIdx = idx

# print("seq max: {} {} {}".format(maxCnt, maxSeq, maxSeqIdx))


## Read the imput output file 
n = int(steinerFile.readline().strip())
treeNodes = []
cells = []
maxAc = -1
for i in range(n):
        l = steinerFile.readline().strip()
        x = l.split()
        treeNodes.append(x[0])
        if x[1] == '1':
                cells.append(x[0])
        ac = sum([1 for s in x[2] if s == 'A'])
        # print("  mx: ac={} len={} x[2]={} i={}".format(ac, len(x[2]), x[2], i))
        if ac > maxAc:
                maxAc = ac
                maxAcLen = len(x[2])
                treeRootSeqIdx = int(x[0])

# print("max: ac:{} ac-len:{} idx:{} ".format(maxAc, maxAcLen, treeRootSeqIdx))

# print("excluding {} {} {}".format(maxAcLen == maxAc, 5 < len(sys.argv), sys.argv[5] == "-exclude-root"))
seqIdxToExclude = -1
if maxAcLen == maxAc and 5 < len(sys.argv) and sys.argv[5] == "-exclude-root":
        seqIdxToExclude = treeRootSeqIdx
        # print("excluding {}".format(seqIdxToExclude))
        
cells = set(cells)

# print("maxSeqIdx: {}".format(maxSeqIdx), file=sys.stderr)
# maxSeqClone = -1

for t in treeNodes:
        print(t, end=" ", file=cloneFile)
        if t in cells and t != seqIdxToExclude and int(t) != seqIdxToExclude:
                print(int(t)+1, end="", file=cloneFile)
        # print("t:{} treeRootSeqIdx:{}".format(t, treeRootSeqIdx), file=sys.stderr)
        if int(t) == treeRootSeqIdx:
                treeRootCloneIdx = int(t)
                # print("  treeRootSeqIdx: {} t: {} treeRootCloneIdx: {}".format(treeRootSeqIdx, t, treeRootCloneIdx), file=sys.stderr)
        print(file=cloneFile)
# print("maxSeqClone: {}".format(maxSeqClone), file=sys.stderr)


print(' '.join(treeNodes), file=treeFile)

edges = {}

m = int(steinerFile.readline().strip())
for i in range(m):
        l = steinerFile.readline().strip()
        x = l.split()
        if x[0] not in edges: edges[x[0]] = []
        if x[1] not in edges: edges[x[1]] = []
        edges[x[0]].append((x[1], x[2]))
        edges[x[1]].append((x[0], x[2]))

# print("edges", end = " ")
# print(edges)

mark = {}
dfsNumCounter = 0
dfsNum = {}
def dfs(v):
        global dfsNumCounter, dfsNum
        dfsNum[v] = dfsNumCounter
        dfsNumCounter += 1
        mark[v] = True
        for (u, _) in edges[v]:
                if u not in mark:
                        dfs(u)

# print('maxSeqClone+1: {}'.format(maxSeqClone+1), file=sys.stderr)
# dfs(str(maxSeqClone+1))
dfs(str(treeRootCloneIdx))
        
treeRootCloneIdx
for v, nei in edges.items():
        for (u,w) in nei:
                if v < u:
                        x,y=u,v
                        if dfsNum[v] > dfsNum[u]:
                                x,y = v,u
                        print("{}->{} {}".format(x, y, w), file=treeFile)

print("seqroot: {} cloneroot: {} exclude: {}".format(treeRootSeqIdx, treeRootCloneIdx, seqIdxToExclude), file=sys.stderr)