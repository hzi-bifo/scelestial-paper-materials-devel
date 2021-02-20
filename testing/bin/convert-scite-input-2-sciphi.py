import sys, itertools, re

sciteFile = open(sys.argv[1])
sciphiFile = open(sys.argv[2], "w+")
sciphiCellNamesFile = open(sys.argv[3], "w+")
sciphiCellRepeat = int(sys.argv[4])

seq = []
for l in sciteFile:
        x = l.strip().split()
        for i, v in enumerate(x):
                while i >= len(seq):
                        seq.append([])
                seq[i].append(int(v))


if len(seq) > 0:
	locCnt = len(seq[0])
	for i in range(locCnt):
		line = '{}\t{}\t{}'.format('seq1', (i+1)*100, 'A')
		r = q = ''
		for s in seq:
			if s[i] == 0:
				r = '.'
			elif s[i] == 1:
				r = 'T'
			elif s[i] == 3:
				r = 'N'
			q = '<'
			line = '{}\t{}\t{}\t{}'.format(line, sciphiCellRepeat, r * sciphiCellRepeat, q * sciphiCellRepeat)
		print(line, file=sciphiFile)

for i, s in enumerate(seq):
	print('{}\t{}'.format(i+1, 'CT'), file=sciphiCellNamesFile)
