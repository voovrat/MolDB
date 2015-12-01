function water_PC(HNC,VUA)

g_dat = load('waterRDFs.txt');
hk = d3fft2(g_dat(:,2:end)-1,g_dat(:,1));

[DU,EU]=unit_names;
rhoA = 0.03336;
rhoB = rhoA * unit2unit(DU,'Angstr','Bohr')^-3;

ck = hk./(1+rhoB*hk);

ck0A = ck(1,1) * unit2unit(DU,'Bohr','Angstr')^3;

kT = 300 * unit2unit(EU,'K','kcal/mol');


P= rhoA*kT - 0.5*rhoA^2*kT*ck0A;
PKBar = P * unit2unit(EU,'kcal/mol','J')*1e30/1e5/1e3;
PC = 0.5*kT*rhoA^2*ck0A*VUA;
Corr = HNC + PC;

fprintf('ck(0): %f [A^3] \nP: %f [kcal*mol^-1/A^3] = %g [KBar]\nPC+: %f [kcal/mol]\nHNC/PC+: %f [kcal/mol]\n',ck0A,P,PKBar,PC,Corr);
