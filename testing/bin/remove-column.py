import sys

toBeRemovedColumns = sys.argv[1].split(',')

first=True
for l in sys.stdin:
	l = l[:-1]
	if first: names = l.split('\t')
	lo=[]
	for n, v in zip(names, l.split('\t')):
		if n not in toBeRemovedColumns: lo.append(v)
	print('\t'.join(lo))
	first=False
