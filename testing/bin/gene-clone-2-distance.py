import sys, itertools

genotypeFile = open(sys.argv[1])
cloneFile = open(sys.argv[2])

fl = genotypeFile.readline()
names = fl.strip().split(' ')
n = len(names)
seq = {name:[] for name in names}
# print(fl)
for l in genotypeFile:
    x = l.strip().split(' ')
    # print(l)
    for (i, v) in enumerate(x[1:]):
        seq[names[i]].append(v)
# print(seq[1])

dist = { i:{ j:-1 for j in names} for i in names }
for (i, s) in seq.items():
    for (j, q) in seq.items():
        x = sum(x != y for x, y in itertools.zip_longest(s, q))
        dist[i][j] = x
        # print(x, end=" ")
#     print()

# print(dist)

clone = {}
for l in cloneFile:
        x = l.strip().split(' ')
        if len(x) == 0: continue
        for se in x[1:]:
                clone[se] = x[0]

# print(clone)

print('""', end=" ")
for (j,d) in clone.items():
        print('"{}"'.format(j), end=" ")
print()
for (i,c) in clone.items():
        print('"{}"'.format(i), end=" ")
        for (j,d) in clone.items():
                print(dist[c][d], end=" ")
        print()
