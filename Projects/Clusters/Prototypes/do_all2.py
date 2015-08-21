
from lmpdriver import *
import sys
import math
import os


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


sc_p = bulk_p.nearliest_neighb(pi,8);

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




sc_rest_oxylist = sc_p.oxylist()[len(c_p.oxylist()):]
sc_rest = sc_p.select_oxygens(sc_rest_oxylist)


def energy_list_rest(c_p,cn,sc_rest,L):
   dEn = []
   for ox in sc_rest.oxylist():
      w = sc_rest.select_oxygens([ox])
      En = L.run(w + cn)
      Ep = L.run(w + c_p)
      dEn.append(Ep-En)
   return dEn


#               []
#            []    [1] 
#       []   [2] [1]    [12]
#    []   [3] [2] [1]  [23][13][12]  [123]  
# 
#   combinations(m,n) =  combinations(m,n-1) + 
#   + combinations(m-1,n-1).append(n)
#
def combinations(m,n):
   if m==0:
      return [[]];
   if m==n: 
      return [range(n)];
   C = combinations(m,n-1)
   C1 = combinations(m-1,n-1)
   for c in C1:
      c.append(n-1)
   return C+C1

def Cmn(m,n):
   if m==0:
      return 1
   if m==n:
      return 1
   return Cmn(m-1,n-1) + Cmn(m,n-1)

#              1
#            1  1
#          1  2  1
#        1  3  3  1
#      1  4  6  4  1 
#    1  5  10 10  5  1
#  1  6  15 20 15  6  1


def energy_list_rest_many(c_p,cn,sc_rest,L,nn):
   dEnn = []
   N = len(sc_rest.oxylist());
   oxl=sc_rest.oxylist()
   cc = combinations(nn,N)
   print 'combinations:'+str(len(cc))
   oxcc = []
   for cox in cc:
      ocox = []
      for j in cox:
         ocox.append(oxl[j])
      print ocox
      oxcc.append(ocox)
      ww = sc_rest.select_oxygens(ocox)
      #print ww
      En = L.run(ww + cn)
      Ep = L.run(ww + c_p)
      dEnn.append(Ep-En)
   return oxcc,dEnn


def nansum(l):
   S = 0
   for x in l:
      if not math.isnan(x):
         S += x
   return S

def nanmean(l):
   S=0
   N=0
   for x in l:
      if not math.isnan(x):
         S+=x
         N+=1
   if N>0:
      return S/N
   else:
      return float('nan');



dEnn = []
dEn = []


for nn in range(1,8):
   print "Rest Cluster Size: "+str(nn)
   oxcc,dE_list = energy_list_rest_many(c_p,c1,sc_rest,L,nn)
   dEnn.append(dE_list)
   dEn.append( sum(dE_list) )


##########

f=open('dEnn.dat');
S1 = f.readline();
f.close();

dEnn = eval(S1);



E_cn = L.run(c1)
E_cp = L.run(c_p)

dE=[]

for nn in range(1,8):
   de=(sum(dEnn[nn-1]) - len(dEnn[nn-1])*(E_cp-E_cn))/Cmn(nn-1,6)
   dE.append(de)


def save_list(l,fname):
   f=open(fname,'w')
   for x in l:
      f.write(str(x)+'\n')
   f.close()


save_list(dE,'dE.dat')

def octave_cmd(s):
   f=open('cmd.m','w')
   f.write(s)
   f.close()


S="""
dE=load('dE.dat');
hnd=plot(dE)
xlabel('neutral cluster size')
ylabel('\Delta E_{p} [eV]')
set(hnd,'LineWidth',4)
print -F:18 dEc.png
"""

octave_cmd(S)

os.system('eog dEc.png')


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

def coul_list_rest_many(c_p,cn,sc_rest,L,nn):
   dEnn = []
   N = len(sc_rest.oxylist());
   oxl=sc_rest.oxylist()
   cc = combinations(nn,N)
   print 'combinations:'+str(len(cc))
   oxcc = []
   for cox in cc:
      ocox = []
      for j in cox:
         ocox.append(oxl[j])
      print ocox
      oxcc.append(ocox)
      ww = sc_rest.select_oxygens(ocox)
      #print ww
      spce_ff(ww)
      hydronium_ff(c_p,-0.5)
      spce_ff(cn)
      En = coulomb(ww,cn)
      Ep = coulomb(ww,c_p)
      dEnn.append(Ep-En)
   return oxcc,dEnn

def addsc(ls,x):
   r = []
   for i in range(len(ls)):
      r.append(ls[i] + x)
   return r

def subsc(ls,x):
   r = []
   for i in range(len(ls)):
      r.append(ls[i] - x)
   return r




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

#dV11/dx_at =  Fbulk + F_c - F_c1
#dC11/dx_p  = F_c 

def energy_rest(dEn,E_c,E_cn):
   dE = 0
   for dEi in dEn:
      dE += dEi - E_c + E_cn
   return dE

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


-2.4 2.2998043395 2.23979909813 
2.2998043395 -2.38 2.26920690991 
2.23979909813 2.26920690991 -2.49 

V =

  -2.4000   2.2998   2.2398
   2.2998  -2.3800   2.2692
   2.2398   2.2692  -2.4900

octave:2> eig(V)
ans =

  -4.7052
  -4.6814
   2.1166

  -2.4000   2.1360   2.2347
   2.1360  -2.7000   2.1029
   2.2347   2.1029  -2.5000

[x,L]=eig(V)
x =

  -0.099437  -0.809000  -0.579337
   0.743604   0.326458  -0.583505
  -0.661184   0.488820  -0.569114


def dd(a,i,j):
   dx = a.x[i] - a.x[j]
   dy = a.y[i] - a.y[j]
   dz = a.z[i] - a.z[j]
   return sqrt( dx*dx + dy*dy + dz*dz )


dd(ideal_hydr,0,1)   # 1.114408046159036
dd(ideal_hydr,0,2)   # 1.0213931788694297
dd(ideal_hydr,0,3)   # 1.0108034948104405

dOH1=sqrt(ideal_hydr.x[1]-ideal_hydr.x[0])**2 + \
ideal_hydr.y[1]-ideal_hydr.y[0])**2 + \
ideal_hydr.z[1]-ideal_hydr.z[0])**2)

