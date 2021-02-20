from PIL import Image
import numpy as np
import sys
im = Image.open(sys.argv[1])

def col(p, i, j):
        gray = [170, 170, 170]
        white = [255, 255, 255]
        blue = [22, 22, 178]
        colors = {'G': gray, 'W': white, 'B': blue}
        #dist = {n: sum([abs(c[k] - p[i,j,k]) for k in range(p.shape[2])]) for n, c in colors.items()}
        if p[i,j,2] > 2 * p[i,j,0] and p[i,j,2] > 2 * p[i,j,0]:
                r = 'B'
        elif all([p[i,j,k] > 200 for k in range(3)]):
                r = 'W'
        else:
                r = 'G'
        return r, p[i,j]

        #r, g, b = p[i,j,0], p[i,j,1], p[i,j,2]
        #isw = all([230 <= p[i,j,k] <=255 for k in range(p.shape[2])])  
        #if isw: return 'www/www/www'
        #return '/'.join([str(p[i,j,k]) for k in range(p.shape[2])])    

colo, cold = [], []
p = np.array(im)
for i in range(p.shape[0]):
        colo.append([])
        cold.append([])
        for j in range(p.shape[1]):
                c, d = col(p, i, j)
                cold[i].append(d)
                colo[i].append(c)
                #print(col(p, i, j), end=" ")
        #print()
for i in range(p.shape[0]):
        for j in range(2, p.shape[1]-2):
                if colo[i][j] == 'G' and (colo[i][j+1] == 'B' or colo[i][j+2] == 'B') and (colo[i][j-1] == 'B' or colo[i][j-2] == 'B'):
                        colo[i][j] = 'W'

uniq_column = []
imp_col_all = []
for j in range(p.shape[1]):
        colu = 0
        for i in range(p.shape[0]):
                if colo[i][j] == 'W':
                        colu += 1
        if colu / p.shape[0] > 0.5:
                colu = 'W'
        else:
                colu = '*'
        if colu == '*':
                imp_col_all.append(j)
        uniq_column.append(colu)

def all_2_range(a):
	imp_col = []
	last_imp_col = -100
	for i in a:
		if i > last_imp_col + 1:
			imp_col.append([i, i+1])
		else:
			imp_col[-1][1] = i+1
		last_imp_col = i
	return imp_col

imp_col = all_2_range(imp_col_all)
	

uniq_row = []
imp_row_all = []
for i in range(p.shape[0]):
        colu = 0
        for j in range(p.shape[1]):
                if colo[i][j] == 'W':
                        colu += 1
        if colu / p.shape[1] > 0.8:
                colu = 'W'
        else:
                colu = '*'
        if colu == '*':
                imp_row_all.append(i)
        uniq_row.append(colu)


imp_row = all_2_range(imp_row_all)

tab = []
for row in imp_row:
	tab.append([])
	for col in imp_col:
		l = []
		for i in range(row[0], row[1]):
			for j in range(col[0], col[1]):
				l.append(colo[i][j])
		c = max(set(l), key = l.count)
		#print('{}:{} {} {}'.format(row, col, c, l))
		tab[-1].append(c)

def col2scite(c):
	if c == 'W': return 3
	if c == 'B': return 1
	return 0

for r in tab:
	for c in r:
		print(col2scite(c), end=' ')
	print()

#print(len(imp_col))
#print(len(imp_row))

#print('{}'.format(' '.join(uniq_column)))
#print('{} {}'.format('-', '/???,???,??? '.join(uniq_column)))
#for i in range(p.shape[0]):
#        print(uniq_row[i], end=' ')
#        for j in range(p.shape[1]):
#                print('{}/{}'.format(colo[i][j], ','.join(["{:03d}".format(cold[i][j][k]) for k in range(3)])), end=' ')
#       # print()
