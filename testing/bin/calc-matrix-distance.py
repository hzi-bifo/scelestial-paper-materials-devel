import sys, os

def loadDist(fName):
	if not os.path.exists(fName):
		return None, None
	#print('loading {}'.format(fName), file=sys.stderr)
	f = open(fName)
	names = f.readline().strip().split()[1:]
	if len(names) == 0:
		return None, None
#	 n = len(names)
	dist = {i:{j:-1 for j in names} for i in names}

	for l in f:
		x = l.strip().split()
		for (i,v) in enumerate(x[1:]):
			dist[x[0]][names[i]] = float(v)
	
	return dist, names

def normalize(d, name):
	# print(d)
	# for r in d:
	#	 print(r)
	s = sum([sum(r.values()) for (i,r) in d.items()])
#	 print(d)
#	 print("sum: {} name: {}".format(s, name))
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

normalizeMethod = 1
argc = 1
if sys.argv[argc] == '-no-normal':
	argc += 1
	normalizeMethod = 0

for i in range(argc, len(sys.argv), 2):
	fn = sys.argv[i]
	d,n = loadDist(fn)
	D.append(d)
	N.append(n)
	T.append(sys.argv[i+1])


NnotNone = [i for i in range(len(N)) if N[i] is not None]
for x in range(len(NnotNone)-1):
	i, j = NnotNone[x], NnotNone[x+1]
	if not set(N[i]) == set(N[j]):
		print("sets {} and {} are non-equal ".format(T[i], T[j]) )
		print(N[i])
		print(N[j])
		raise("invalid set of names {}, {}".format(N[i], N[j]))

DNormal = []
for (i,d) in enumerate(D):
	if d is not None:
		DNormal.append(normalize(d, T[i]))
	else:
		DNormal.append(None)

	# print("{}".format(T[i]))
	# for (n,l) in DNormal[-1].items():
	#		 print("n: {}".format(n), end=" ")
	#		 print(l)

if normalizeMethod == 1:
	D = DNormal

for j in range(len(T)):
	print(T[j], end=" ")
print()
	
for i in range(len(T)):
	print(T[i], end=" ")
	for j in range(len(T)):
		if N[i] is not None and N[j] is not None:
			d = "{0:.4f}".format(distance(D[i], D[j], N[NnotNone[0]]))
		else:
			d = 'NA'
		print(d, end=" ")
	print()

