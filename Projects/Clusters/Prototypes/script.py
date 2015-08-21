from lmpdriver import *
import sys

bulk_p = XYZ()
bulk_p.read('oktamer_proton.xyz');

pi = bulk_p.proton_index();

bulk_p2 = bulk_p.nearliest_neighb(pi,8);

pi1 = bulk_p2.proton_index()
pi2 = 1
pi3 = 2

bulk1 = bulk_p2.copy()
bulk1.remove(pi1)

bulk2 = bulk_p2.copy()
bulk2.remove(pi2)

bulk3 = bulk_p2.copy()
bulk3.remove(pi3)

c_p = bulk_p.nearliest_neighb(pi,1);
c1 = c_p.copy()
c1.remove(c_p.proton_index())

c2 = c_p.copy()
c2.remove(1)

c3 = c_p.copy()
c3.remove(2)


sc_p = bulk_p.nearliest_neighb(pi,4);

sc1 = sc_p.copy()
sc1.remove(sc_p.proton_index())

sc2 = sc_p.copy()
sc2.remove(1)

sc3 = sc_p.copy()
sc3.remove(2)

ideal_hydr = bulk_p.nearliest_neighb(pi,1)
ideal_h2o = ideal_hydr.copy()
ideal_h2o.remove(ideal_hydr.proton_index())

L = LMP('../../src/lmp_debug')

E_h2o = L.run(ideal_h2o)
E_hydr = L.run(ideal_hydr)
E_p0 = E_hydr - E_h2o

E_bulk1 = L.run(bulk1) 
E_bulk2 = L.run(bulk2) 
E_bulk3 = L.run(bulk3) 


E_c = L.run(c_p)
E_c1 = L.run(c1)
E_c2 = L.run(c2)
E_c3 = L.run(c3)


V11 = E_bulk1 + (E_c - E_c1) #- 8*E_h2o - E_p0
V22 = E_bulk2 + (E_c - E_c2) #- 8*E_h2o - E_p0
V33 = E_bulk3 + (E_c - E_c3) #- 8*E_h2o - E_p0

E_sc = L.run(sc_p)
E_sc1 = L.run(sc1)
E_sc2 = L.run(sc2)
E_sc3 = L.run(sc3)

dE1 = E_sc - E_sc1 - (E_c - E_c1)
dE2 = E_sc - E_sc2 - (E_c - E_c2)
dE3 = E_sc - E_sc3 - (E_c - E_c3)

V12 = sqrt(dE1*dE2)
V13 = sqrt(dE1*dE3)
V23 = sqrt(dE2*dE3)

V=[ [ V11, V12, V13],
  [ V12, V22, V23],
  [ V13, V23, V33 ] ]

for i in range(3):
   for j in range(3):
      sys.stdout.write(str(V[i][j]) + " ")
   sys.stdout.write("\n")

