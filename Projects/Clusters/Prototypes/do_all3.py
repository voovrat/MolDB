from lmpdriver import *

c8p=XYZ()
c8p.read('oktamer_proton.xyz')
c8n=c8p.unselect([c8p.proton_index()]);

def spce_ff(cwater):
   cwater.charge = []
   for i in range(len(cwater)):
      if cwater.atoms[i] == 'O':
         ch = -0.8476
      else:
         ch = 0.4238
      cwater.charge.append(ch)

def hydronium_ff(c,qH):
   c.charge = [ 1 - 3*qH, qH, qH, qH]


def coulomb(c1,c2):
   U=0
   for i in range(len(c1)):
      x1 = c1.x[i]
      y1 = c1.y[i]
      z1 = c1.z[i]
      q1 = c1.charge[i]
      for j in range(len(c2)):
         B = 0.52918 #  Angstr -> Bohr
         dx = (c2.x[j] - x1)/B
         dy = (c2.y[j] - y1)/B
         dz = (c2.z[j] - z1)/B
         q2 = c2.charge[j]
         r = sqrt( dx*dx + dy*dy + dz*dz)
         U += q1*q2/r
   return U*27.211385   # Hartree -> eV




def qH_best(c8p,nc):
   c = c8p.nearliest_neighb(c8p.proton_index(),1)
   cc = c8p.nearliest_neighb(c8p.proton_index(),nc);
   c_rest_oxy = cc.oxylist()[1:]
   c_rest = cc.select_oxygens(c_rest_oxy)
   pass
   spce_ff(c_rest)
   hydronium_ff(c,0.5)
   c1n = c.select_oxygens([1])
   pass
   pass
   L=LMP('../../src/lmp_debug')
   pass
   E_ccp = L.run(c_rest + c)
   pass
   E_c1p = L.run(c)
   E_c1n = L.run(c1n)
   E_rest = L.run(c_rest)
   pass
   U_real = E_ccp - E_rest - E_c1p
   U_bond = E_c1p - E_c1n
   U_coul = coulomb(c,c_rest)
   pass
   qH_list = []
   U_coul_list = []
   dU_list = []
   pass   
   qH = -1
   Nq = 100
   minDelta = 100
   for i in range(2*Nq):
      hydronium_ff(c,qH)
      U_coul = coulomb(c,c_rest)
      dU = U_coul - U_real
      qH_list.append(qH)
      U_coul_list.append(U_coul)
      dU_list.append(dU)
      qH = qH + 1.0/Nq
      if abs(dU) < minDelta:
         best_q = qH
         minDelta = abs(dU)
   pass
   return best_q

for i in range(2*Nq):
   print qH_list[i],dU_list[i]


q_best_list = []

for i in range(2,9):
   q_best = qH_best(c8p,i)
   q_best_list.append(q_best)


def save_list(l,fname):
   f=open(fname,'w')
   for x in l:
      f.write(str(x)+'\n')
   f.close()

def octave_cmd(s):
   f=open('cmd.m','w')
   f.write(s)
   f.close()


save_list(q_best_list,'q_best_list.dat')

s="""
q_best_list=load('q_best_list.dat')
plot(2:8,q_best_list)
"""
octave_cmd(s)




cp  rest

c sc bulk

c + bulk

bulkp \approx bulk_n + bond + dcoul


dE = (cp -- rest) - (cn -- rest) 

cp -- rest = sum(w in rest) (cp+w) - (cn+w)

p+H2O

bond = cp - cn
dcoul = coul(cp,bulk_rest) - coul(cn,bulk_rest)

sc 


########

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


%%%%%%%%%%%%%%%%%%%%


E_cp = load('E_cp.dat')
E_cn = load('E_cn.dat')
E_coul = load('E_coul.dat')

dE = E_cp - E_cn - E_coul
dE0 = E_cp - E_cn

hnd=plot([dE-dE(end) dE0-dE0(end)])
legend('E_c^p - E_c^n - E_{coul}','E_c^p - E_c^n')

set(hnd,'LineWidth',4)
xlabel('Cluster Size')
ylabel('\Delta E [ eV]')

print dE.png -F:18

[DU,EU]=unit_names

hnd=plot((1:8),[dE-dE(end) dE0-dE0(end)]*unit2unit(EU,'eV','kcal/mol'),[1 8],[5 5],[1 8],[-5 -5])
legend('E_c^p - E_c^n - E_{coul}','E_c^p - E_c^n','5 kcal/mol interval')

set(hnd,'LineWidth',4)
set(hnd([3 4]),'LineWidth',2,'Color',[0 0 0],'LineStyle','--')
xlabel('Cluster Size')
ylabel('\Delta E [ kcal/mol]')

print dE_kcal_mol.png -F:18







