import sys
from lmpdriver import *

if len(sys.argv)<5:
   print  'Usage: python gen_cluster.py  c.xyz sc_rest.xyz subindeces(in format id1%id2%..%idN )  output.xyz'
   quit()

c=XYZ()
c.read(sys.argv[1])

sc_rest = XYZ()
sc_rest.read(sys.argv[2])

words = str.replace(sys.argv[3],'c',' ').split()
oxl = [] #sc_rest.oxylist()
for w in words:
   oxl.append( sc_rest.oxylist()[int(w)])

cc = c + sc_rest.select_oxygens(oxl)
cc.write(sys.argv[4])

