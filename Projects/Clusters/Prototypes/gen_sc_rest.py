from lmpdriver import *
a = XYZ()
a.read('oktamer_proton.xyz')
c8p=a.nearliest_neighb(a.proton_index(),8)
sc_rest = c8p.select_oxygens( c8p.oxylist()[1:])
sc_rest.write('sc_rest.xyz')
