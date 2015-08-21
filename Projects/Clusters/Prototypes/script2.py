
from lmpdriver import *

c8p=XYZ()
c8p.read('oktamer_proton.xyz')
c8n=c8p.unselect([c8p.proton_index()]);

cp=[]
cn=[]

for i in range(8):
   cpi=c8p.nearliest_neighb(c8p.proton_index(),i+1)
   cni = cpi.unselect([cpi.proton_index()]);
   cp.append(cpi)
   cn.append(cni)

E_cp = []
E_cn = []
E_coul = []

L = LMP('../../src/lmp_debug')

for i in range(8):
   e_cp = L.run(cp[i])
   e_cn = L.run(cn[i])
   e_coul = cp[i].proton_coulomb(cp[i].proton_index())
   E_cp.append(e_cp)
   E_cn.append(e_cn)
   E_coul.append(e_coul)

def save_list(l,fname):
   f=open(fname,'w')
   for x in l:
      f.write(str(x)+'\n')
   f.close()

save_list(E_cp,'E_cp.dat')
save_list(E_cn,'E_cn.dat')
save_list(E_coul,'E_coul.dat')


