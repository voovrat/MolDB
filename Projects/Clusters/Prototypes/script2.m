
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




