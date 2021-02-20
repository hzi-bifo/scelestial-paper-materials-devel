import sys, re, argparse
from graphviz import Digraph
from lib import *

parser = argparse.ArgumentParser(description='Clone tree to mutation tree.')
parser.add_argument('inputFileName')
parser.add_argument('outFileName')
#parser.add_argument('--has-row-head', dest='hasRowHead', action='store_const', const = False, default = True)
#parser.add_argument('--has-col-head', dest='hasColHead', action='store_const', const = False, default = True)
parser.add_argument('--row', dest='rowMVThreshold', action='store', type=float, default = 0)
parser.add_argument('--col', dest='colMVThreshold', action='store', type=float, default = 0)

args = parser.parse_args()

sequences = loadSequenceFile(args.inputFileName)

newSeq = [ seq for seq in sequences if sum([x != 3 for x in seq]) >= len(seq) * args.colMVThreshold ]
sequences = newSeq
#print(sequences)
print("{}".format(", ".join([str(len(seq)) for seq in sequences])))
if len(sequences) > 0:
	validColumns = [ i for i in range(len(sequences[0])) if sum([ seq[i] != 3 for seq in sequences]) >= len(sequences) * args.rowMVThreshold]
	validColumns = set(validColumns)
	newSeq = [ [x for i,x in enumerate(seq) if i in validColumns ] for seq in sequences]
	sequences = newSeq
	print("{}".format(", ".join([str(len(seq)) for seq in sequences])))

if len(sequences) > 0:
	print("{} x {}".format(len(sequences[0]), len(sequences)))
writeSequenceFile(sequences, args.outFileName)
