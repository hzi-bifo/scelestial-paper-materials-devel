import sys
from pygraphml import GraphMLParser
from pygraphml import Graph

parser = GraphMLParser()
g = parser.parse(sys.argv[1])

ftree = open(sys.argv[2], "w")
fclone = open(sys.argv[3], "w")

for v in g.nodes():
    l = v['v_members']
    ll = [str(int(x)+1) for x in l.split()]
    l = ' '.join(ll)
    fclone.write("{0} {1}\n".format(v.id, l))

for v in g.nodes():
    ftree.write("{} ".format(v.id))
ftree.write("\n")

#g.show()
# for v in g.nodes():
#     print("{}: ".format(v.id), end=" ")
#     for u in v.children():
#         print(u.id, end = " ")
#     print()

for e in g.edges():
    ftree.write("{}->{} {}\n".format(e.node2.id, e.node1.id, e['e_weight']))
        
