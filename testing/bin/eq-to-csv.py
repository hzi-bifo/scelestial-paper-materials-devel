import sys, re, argparse, math

maxattrcount = 20

parser = argparse.ArgumentParser(description='Generate samples from single-cell data set.')
parser.add_argument('--avg', dest='avg1', action='store', required = False, nargs='*')
parser.add_argument('--var', dest='var1', action='store', required = False, nargs='*')
parser.add_argument('--time', dest='time', action='store', required = False, nargs='*')
for i in range(2, maxattrcount):
	parser.add_argument('--avg'+str(i), dest='avg'+str(i), action='store', required = False, nargs='*')
	parser.add_argument('--var'+str(i), dest='var'+str(i), action='store', required = False, nargs='*')
args = parser.parse_args()

def msq(x):
	return x * x

def mean(v):
	return sum(v) / len(v)

def var(v):
	m = mean(v)
	if len(v) == 1: return 0
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
for i in range(1,maxattrcount):
	x = getattr(args, "avg" + str(i))
	if x is not None:
		merges.append((x, mean))
	x = getattr(args, "var" + str(i))
	if x is not None:
		merges.append((x, var))

for dlist, meth in merges:
	if dlist is not None:
		vars = set(vars) | set([dlist[0]])

vars = sorted(vars)

print("\t".join(vars))
for ll in lines:
	for v in vars:
		for dlist, meth in merges:
			if dlist is not None and v == dlist[0]:
				vl = [float(ll[mg]) for mg in dlist[1:]]
				ll[v] = meth(vl)
		if v in varTime:
			m = re.match("(\d+)m(\d+\.\d+)s", ll[v])
			ll[v] = int(m[1]) * 60 + float(m[2])
		print(ll[v], end="\t")
	print()
	
