import sys
import statistics

args_samples_set = [50, 100]
args_clones_set = [5, 10]
args_sites_set = [20, 50]

def read_int_args(argv, argc):
	r = []
	while argc < len(argv) and argv[argc][0] != '-':
		r.append(int(argv[argc]))
		argc += 1
	return argc, r

calculate_variance = False
argc = 1
while argc < len(sys.argv):
	if sys.argv[argc][0] == '-':
		if sys.argv[argc] == '--samples':
			argc, args_samples_set = read_int_args(sys.argv, argc+1)
		if sys.argv[argc] == '--clones':
			argc, args_clones_set = read_int_args(sys.argv, argc+1)
		if sys.argv[argc] == '--sites':
			argc, args_sites_set = read_int_args(sys.argv, argc+1)
		if sys.argv[argc] == '--var':
			calculate_variance = True
			argc += 1

def isnumeric(s):
	'''returns True if string s is numeric'''
	return all(c in "0123456789.+-" for c in s) and any(c in "0123456789" for c in s)

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
	for samples in args_samples_set:
		for clones in args_clones_set:
			for sites in args_sites_set:
				elements = []
				for val in values:
					if val['sample'] == str(samples) and val['clone'] == str(clones) and val['site'] == str(sites):
						elements.append(float(val[alg])) # +':'+str(samples)+str(clones)+str(sites))
				l.append('{:.2f}'.format(sum(elements) / len(elements))) # + ' (' + str(len(elements)) + ')')
				if calculate_variance:
					l[-1] += ' ({:.2f})'.format(statistics.variance(elements))
	lines.append(l)

def first_number(val):
	if val.find('(') != -1: val = val[:val.find('(')]
	if val.find(' ') != -1: val = val[:val.find(' ')]
	return val

mincol = []
for i in range(len(lines[0])):
	mini = '100'
	for l in lines:
		val = first_number(l[i])
		#print(val, file=sys.stderr)
		if l[i] != 'nan' and isnumeric(val) and float(val) < float(mini):
			mini = val
	mincol.append(mini)
print(mincol)
for j, alg in enumerate('onco,impt,bphl,scte,sasc,sciphi,sifit,siclonefit'.split(',')):
	print('&'.join([alg] + [x if first_number(x) != mincol[i] else '{\\bf ' + x + '}' for i, x in enumerate(lines[j])]) + '\\'*2)

