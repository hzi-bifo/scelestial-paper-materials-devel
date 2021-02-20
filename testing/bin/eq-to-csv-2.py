import sys, re, argparse, math

parser = argparse.ArgumentParser(description='Generate samples from single-cell data set.')
parser.add_argument('-a', dest='arg', action='append', required = False, nargs=3)
parser.add_argument('--len', dest='length', action='store', required = True, type=int)
parser.add_argument('--time', dest='time', action='store', required = False, nargs='*')
parser.add_argument('-by-alg', '--byAlg', action='store_true')
parser.add_argument('-mean-na', dest='meanNA', default=0)
args = parser.parse_args()

def msq(x):
	return x * x

def mean(v):
	if len(v) == 0: return args.meanNA
	return sum(v) / len(v)

def var(v):
	m = mean(v)
	if len(v) <= 1: return 0
	return math.sqrt(sum([msq(x-m) for x in v]) / (len(v)-1.0) ) 

vars = []
lines = []
for l in sys.stdin:
	x = re.split(r'\s', l.strip())
	ll = {}
	for xx in x:
		#print(xx)
		vr, vl = xx.split(':')
		ll[vr] = vl
	lines.append(ll)
	vars = set(vars) | set(ll.keys())

varTime = set([])
if args.time is not None:
	varTime = set(args.time)
#merges = [(args.avg1, mean), (args.avg2, mean), (args.avg3, mean), (args.avg4, mean), (args.avg5, mean), (args.avg6, mean), (args.avg7, mean), (args.avg8, mean), (args.var1, var), (args.var2, var), (args.var3, var), (args.var4, var), (args.var5, var), (args.var6, var), (args.var7, var), (args.var9, var)]
merges = []
for a in args.arg:
	dlistAvg = [a[1]]
	dlistVar = [a[2]]
	for i in range(args.length):
		s = a[0].replace('?', str(i))
		dlistAvg.append(s)
		dlistVar.append(s)
	merges.append((dlistAvg, mean))
	merges.append((dlistVar, var))

for dlist, meth in merges:
	if dlist is not None:
		vars = set(vars) | set([dlist[0]])

vars = sorted(vars)

def isfloat(s):
	try:
		float(s)
		return True
	except:
		return False

if args.byAlg: 
	for ll in lines:
		for a in args.arg:
			line = [a[0]]
			for i in range(args.length):
				s = a[0].replace('?', str(i))
				line.append(ll[s] if s in ll else "NA")
			print('\t'.join(line))
else:
	print("\t".join(vars))
	for ll in lines:
		for v in vars:
			for dlist, meth in merges:
				if dlist is not None and v == dlist[0]:
					vl = [float(ll[mg]) if mg in ll and isfloat(ll[mg]) else "NA" for mg in dlist[1:]]
					ll[v] = meth([v for v in vl if v != "NA"])
			if v in varTime:
				m = re.match("(\d+)m(\d+\.\d+)s", ll[v])
				ll[v] = int(m[1]) * 60 + float(m[2])
			print(ll[v], end="\t")
		print()
	
