import sys, re, argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--rep', dest='rep', action='store', required = True)
parser.add_argument('--pat1', dest='pat1', action='store', required = False, nargs='*')
parser.add_argument('--pg1', dest='pg1', action='store', required = False, nargs='*')
parser.add_argument('--pg2', dest='pg2', action='store', required = False, nargs='*')
parser.add_argument('--sg1', dest='sg1', action='store', required = False, nargs='*')
parser.add_argument('--sg1-fault', dest='sg1Fault', action='store', required = False, nargs=1)
parser.add_argument('--print', dest='print', action='store', required = False, nargs='*')
args = parser.parse_args()

def tosecond(s):
	m = re.search("([0-9]+)m([0-9.]+)s", s)
	if m:
		return int(m.group(1)) * 60 + float(m.group(2))
	return float(s)
#patterns = [ ["^DIR=", "SAMPLE=", "sample", "LOCUS=", "locus", "STEP=", "step", "MVR=", "mvr", "0->1=", "zor", "1->0:", "ozr", "STEPMUT:", "stepmut", "CLONE:", "clone"] ]
patterns = [ args.pat1 ]
#patGen = [ [args.rep, "True", ["\\s O?", "\\s I?", "\\s B?", "\\s S?"]], [args.rep, "pTrue", ["\\s pO?", "\\s pI?", "\\s pB?", "\\s pS?"]] ]
patGen = [ [args.rep, args.pg1[0], args.pg1[1:]], [args.rep, args.pg2[0], args.pg2[1:]] ]
varSeqs = []
varSeqsFault = []

for i, sg in enumerate([args.sg1]):
	if sg is None: continue
#print("sg: {} {}".format(sg, args.rep))
	r = [sg[0]]
	for i in range(int(args.rep)):
		for s in sg[1:]:
			r.append(s.replace("?", str(i)))
	varSeqs.append(r)
	if args.sg1Fault and i <= 0:
		varSeqsFault.append(args.sg1Fault[0])
	else:
		varSeqsFault.append('a^')

for pg in patGen:
	for i in range(int(pg[0])):
		pat1 = [pg[1].replace("?", str(i))]
		for s in pg[2]:
			for ss in s.replace("?", str(i)).split():
				pat1.append(ss)
		patterns.append(pat1)
	
print("patterns: {}".format(patterns), file=sys.stderr)
print("patterns: {}".format(varSeqs), file=sys.stderr)
varVal = {}
varSeqSeenCount = [0 for vs in varSeqs]
varSeqOcc = [-1 for vs in varSeqs]
varSeqFaultOcc = [-1 for vs in varSeqs]

lineNo = -1
for l in sys.stdin:
	lineNo += 1
	l = l[:-1]
	for printPattern in args.print:
		if re.search(printPattern, l) and len(varVal) > 0:
			for vr, vl in varVal.items():
				print("{}:{}".format(vr, vl), end=" ")
			print()
			varSeqSeenCount = [0 for vs in varSeqs]
			varSeqOcc = [-1 for vs in varSeqs]
			varSeqFaultOcc = [-1 for vs in varSeqs]
			varVal = {}
	for varLine in patterns:
		if re.search(varLine[0], l):
			pat = ""
			for i in range(1, len(varLine), 2):
				pat = pat + " *" + varLine[i] + ' *(\S+)'
			#print("pat: {}".format(pat), file=sys.stderr)
			m = re.search(pat, l)
			if m:
				for i in range(1, len(varLine), 2):
					varVal[varLine[i+1]] = m.group(int(i/2)+1)
			else:
				print("! INVALID PATTERN MATCH !!{}, '{}'".format(l, pat), file=sys.stderr)
	for i, varSeq in enumerate(varSeqs):
		try:
			if re.search(varSeq[0], l):
				if '(' not in varSeq[0]:
					m = re.search(varSeq[0] + '[ \t]*(\S+)', l)
					varVal[varSeq[varSeqSeenCount[i]+1]] = tosecond(m.group(1))
				else:
					m = re.search(varSeq[0], l)
					g = 1
					if len(m.group(1)) == 0: g += 1
					varVal[varSeq[varSeqSeenCount[i]+1]] = tosecond(m.group(g))
				#print('load {} = {} // {}'.format(varSeq[varSeqSeenCount[i]+1], varVal[varSeq[varSeqSeenCount[i]+1]], l), file=sys.stderr)
				
				if not args.sg1Fault:
					varSeqSeenCount[i] += 1
				varSeqOcc[i] = lineNo
			if re.search(varSeqsFault[i], l):
				#print(' F line: {}/{}/ match:{} fmatch:{} IG:{}'.format(lineNo, l, varSeqOcc[i], varSeqFaultOcc[i], varSeqOcc[i] < varSeqFaultOcc[i] and varSeqFaultOcc[i] != -1), file=sys.stderr)
				#if varSeqOcc[i] < varSeqFaultOcc[i]: # and varSeqFaultOcc[i] != -1:
				#	varVal[varSeq[varSeqSeenCount[i]+1]] = 'NA'
				if args.sg1Fault:
					#print('  done {} // {}'.format(varSeq[varSeqSeenCount[i]+1], l), file=sys.stderr)
					if varSeq[varSeqSeenCount[i]+1] not in varVal:
						varVal[varSeq[varSeqSeenCount[i]+1]] = 'NA'
					varSeqSeenCount[i] += 1
				varSeqFaultOcc[i] = lineNo
		except:
			print("Unexpected error:", sys.exc_info()[0], "line: ", l, 'varseq:', varSeq)
			raise
			

