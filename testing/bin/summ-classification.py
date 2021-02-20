import sys

def loadDist(fName):
    f = open(fName)
    names = [int(s.replace('"', '', 2)) for s in f.readline().strip().split()[1:]]
#     n = len(names)
    dist = {i:{j:-1 for j in names} for i in names}

    for l in f:
        x = l.strip().split()
        for (i,v) in enumerate(x[1:]):
            dist[int(x[0].replace('"', '', 2))][names[i]] = float(v)
    
    return dist, names

def normalize(d, name):
    # print(d)
    # for r in d:
    #     print(r)
    s = sum([sum(r.values()) for (i,r) in d.items()])
#     print(d)
#     print("sum: {} name: {}".format(s, name))
    if s == 0:
            return { i:{ j:v for (j,v) in r.items() } for (i,r) in d.items() }
    return { i:{ j:v/s for (j,v) in r.items() } for (i,r) in d.items() }

def distance(d1, d2, n):
    di = 0
    for i in n:
        for j in n:
            di += abs(d1[i][j] - d2[i][j])
    return di


T, D, N = [], [], []

for i in range(2, len(sys.argv), 2):
        fn = sys.argv[i]
        d,n = loadDist(fn)
        D.append(d)
        N.append(n)
        T.append(sys.argv[i+1])


for i in range(len(N)-1):
    if not set(N[i]) == set(N[i+1]):
        print("sets {} and {} are non-equal ".format(T[i], T[i+1]) )
        print(N[i])
        print(N[i+1])
        raise("invalid set of names {}, {}".format(N[i], N[i+1]))

for i in range(len(N)):
    for i, c in enumerate(N[i]):
            if int(c) != i+1:
                raise("invalid name {}, {}, {}".format(T[i], c, i+1))


DNormal = []
for (i,d) in enumerate(D):
        DNormal.append(normalize(d, T[i]))

        # print("{}".format(T[i]))
        # for (n,l) in DNormal[-1].items():
        #         print("n: {}".format(n), end=" ")
        #         print(l)

D = DNormal


cellClass = []
f = open(sys.argv[1])
for l in f:
        cellClass.append(int(l.strip()))
classSet = set(cellClass)
classCells = {cls:[cell for cell in N[0] if cellClass[cell-1] == cls ] for cls in classSet}

intraClassDist = []
interClassDist = []

def intraClassDistNorm(dist):
        return sum(dist)/len(dist)

def interClassDistNorm(dist):
        return min(dist)
        # return sum(dist)/len(dist)

for alg, (dist, names, title) in enumerate(zip(D, N, T)):
        intraClassDist.append({})
        interClassDist.append({})
        for cls in classSet:
                dist = []
                for cell1 in classCells[cls]:
                        for cell2 in classCells[cls]:
                                dist.append(D[alg][cell1][cell2])
                intraClassDist[-1][cls] = intraClassDistNorm(dist)

        for cls1 in classSet:
                for cls2 in classSet:
                        dist = []
                        for cell1 in classCells[cls1]:
                                for cell2 in classCells[cls2]:
                                        dist.append(D[alg][cell1][cell2])
                        interClassDist[-1][(cls1, cls2)] = interClassDistNorm(dist)

# print(intraClassDist)

print("Intra Class Distances:")
print("{}".format(' '.join([str(cls) for cls in classSet])))
for alg, (dist, names, title, icDist) in enumerate(zip(D, N, T, intraClassDist)):
        print("{} {}".format(title, ' '.join(["{0:.6f}".format(dist) for cls, dist in icDist.items()])))

print("Inter Class Distances:")
for alg, (dist, names, title, icDist) in enumerate(zip(D, N, T, intraClassDist)):
        print(title)
        # print(interClassDist[alg])
        print("{}".format(' '.join([str(cls) for cls in classSet])))
        for cls in classSet:
                print("{} ".format(cls), end="")
                for cls2 in classSet:
                        print("{0:.6f} ".format(interClassDist[alg][(cls, cls2)]), end="")
                print()

# for j in range(len(T)):
#     print(T[j], end=" ")
# print()
    
# for i in range(len(T)):
#     print(T[i], end=" ")
#     for j in range(len(T)):
#         d = distance(D[i], D[j], N[0])
#         print("{0:.2f}".format(d), end=" ")
#     print()

