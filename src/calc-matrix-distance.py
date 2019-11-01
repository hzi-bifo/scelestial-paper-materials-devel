import sys

def loadDist(fName):
	f = open(fName)
	names = f.readline().strip().split()[1:]
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

for i in range(1, len(sys.argv), 2):
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

DNormal = []
for (i,d) in enumerate(D):
		DNormal.append(normalize(d, T[i]))

		# print("{}".format(T[i]))
		# for (n,l) in DNormal[-1].items():
		#		 print("n: {}".format(n), end=" ")
		#		 print(l)


D = DNormal

for j in range(len(T)):
	print(T[j], end=" ")
print()
	
for i in range(len(T)):
	print(T[i], end=" ")
	for j in range(len(T)):
		d = distance(D[i], D[j], N[0])
		print("{0:.2f}".format(d), end=" ")
	print()

