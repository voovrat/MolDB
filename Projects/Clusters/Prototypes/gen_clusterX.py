import sys
from lmpdriver import *

if len(sys.argv)<5:
   print  'Usage: python gen_cluster.py  c.xyz frame.xyz subindeces(in format w<id1>w<id2>w..w<idN> )  output.xyz'
   print 'replaces the 1st molecule by c.xyz, generates the cluster of the rest'
   quit()

c=XYZ()
c.read(sys.argv[1])

fram  = XYZ();
fram.read(sys.argv[2])



words = str.replace(sys.argv[3],'w',' ').split()
oxl = [] #sc_rest.oxylist()
for w in words:
  oxl.append( int(w) )

c0 = fram.select_oxygens([oxl[0]]);
rest  = fram.select_oxygens(oxl[1:]);  

x0 = c0.x[0]
y0 = c0.y[0]
z0 = c0.z[0]

x1 = c.x[0]
y1 = c.y[0]
z1 = c.z[0] 


for i in range(len(c)):
  c.x[i] = c.x[i] - x1 + x0
  c.y[i] = c.y[i] - y1 + y0
  c.z[i] = c.z[i] - z1 + z0

cc = c + rest;

cc.write(sys.argv[4])

