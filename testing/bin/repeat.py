import sys

pat = sys.argv[1]
for i in range(int(sys.argv[2]), int(sys.argv[3])):
	print(pat.replace('?', str(i)), end=" ")
