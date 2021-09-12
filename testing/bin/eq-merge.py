import sys, re, argparse, math

parser = argparse.ArgumentParser(description='Generate samples from single-cell data set.')
parser.add_argument('--on', dest='mergeOn', action='store', required = True, nargs=1)
args = parser.parse_args()
args.mergeOn = args.mergeOn[0]

lines = []
onValues = []
for l in sys.stdin:
	l = l[:-1]
	lv = {x.split(':')[0]:x.split(':')[1] for x in l.split()}
	lines.append(lv)
	onValues.append(lv[args.mergeOn])

#print(lines, file=sys.stderr)
onValues = list(set(onValues))
mergedLines = []
onValueObserved = {ov:0 for ov in onValues}
for lv in lines:
	l = onValueObserved[lv[args.mergeOn]]
	while len(mergedLines) <= l: 
		mergedLines.append({})
	for k, v in lv.items():
		#print('add to merged lines l={} k={} k-replaced={} lml={}'.format(l, k, k.replace('0', lv[args.mergeOn], 1), len(mergedLines)), file=sys.stderr)
		mergedLines[l][k.replace('0', lv[args.mergeOn], 1)] = v
		#print('  new line={}'.format(mergedLines[l]), file=sys.stderr)
	onValueObserved[lv[args.mergeOn]] += 1

for l in mergedLines:
	print(' '.join(k+':'+v for k,v in l.items()))
