import sys
names, values, first = [], [], True
for l in sys.stdin:
	if first:
		names = l.split('\t')
	else:
		values.append(dict(zip(names, l.split('\t'))))
	first = False
lines = []
for alg in 'onco,impt,bphl,scte,sasc,sciphi,sifit,siclonefit'.split(','):
	l = []
	for samples in [50, 100]:
		for clones in [5, 10]:
			for sites in [20, 50]:
				for val in values:
					if val['sample'] == str(samples) and val['clone'] == str(clones) and val['site'] == str(sites):
						l.append(val[alg]) # +':'+str(samples)+str(clones)+str(sites))
	lines.append(l)
mincol = []
for i in range(len(lines[0])):
	mini = '100'
	for l in lines:
		if l[i] != 'nan' and float(l[i]) < float(mini):
			mini = l[i]
	mincol.append(mini)
print(mincol)
for j, alg in enumerate('onco,impt,bphl,scte,sasc,sciphi,sifit,siclonefit'.split(',')):
	print('&'.join([alg] + [x if x != mincol[i] else '{\\bf ' + x + '}' for i, x in enumerate(lines[j])]) + '\\'*2)

