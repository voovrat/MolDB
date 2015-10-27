
from lmpdriver import *
L=LMP('/home/voovrat/460GB/Development/LAMMPS/EVB-NN/lammps/src/lmp_debug')
a = XYZ()
a.read('cluster.xyz')
e=L.run(a)

f=open('energy.dat','w')
f.write( "%0.15f\n" % e )
f.close()

f=open('energy_tot.dat','w')
f.write( "%0.15f\n" % L.etot )
f.close()

