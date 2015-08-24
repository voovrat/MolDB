charges = doFit('Yes','molecules.lst');
nocharges = doFit('No','molecules.lst');
hydro = doFit('No','molecules.lst','MDFT1.chain','','water_model,hydro_type,closure');
%hydro2 = doFit('No','selected.lst','MDFT1.chain',',hs=T,water_model=spc');

plotFit(charges,'WithCharges','MD');
plotFit(nocharges,'NoCharges','MD');
%plotFit(hydro1,'Hydrophobicity','MD');
plotFit(hydro,'Hydrophobicity','MD');

scaling_spc=doFit('No','selected.lst','MDFT2.chain',',water_model=spc','hydro_scaling,closure');
scaling_spce=doFit('No','selected.lst','MDFT2.chain',',water_model=spce','hydro_scaling,closure');

plotFit(scaling_spc,'ChandlerScalingSPC','MD');
plotFit(scaling_spce,'ChandlerScalingSPCE','MD');


vdw_scaling_spc=doFit('No','selected.lst','MDFT_Hydro_VdW.chain',',water_model=spc','a_scaling,closure');
vdw_scaling_spce=doFit('No','selected.lst','MDFT_Hydro_VdW.chain',',water_model=spce','a_scaling,closure');

plotFit(vdw_scaling_spc,'VdWScalingSPC','MD');
plotFit(vdw_scaling_spce,'VdWScalingSPCE','MD');

VdW = doFit('No','molecules.lst','MDFT_Hydro_Both.chain',',hydro_type=VdW,water_model=spc,closure=HNC','water_model,closure,hydro_type');

figure
h=plot(VdW.Exp,VdW.FE,'s',[min(VdW.Exp) max(VdW.Exp)],[min(VdW.Exp) max(VdW.Exp)],'k')
axis equal
set(h,'MarkerSize',3)
set(h,'LineWidth',2);
xlabel('\Delta G_{MD} [kcal/mol]')
ylabel('\Delta G_{calc} [kcal/mol]')

title('VdW/SPC/HNC')

print -dpng VdW.png

[H,X]=hist(VdW.FE - VdW.Exp,100);

dx = X(2)-X(1);

M = nanmean(VdW.FE - VdW.Exp);
S = nanstd(VdW.FE - VdW.Exp);

N = sum(H);

figure
plot(X,H,X,N*dx*gaussian(X,M,S),'Linewidth',3)
xlabel('\Delta G_{calc} - \Delta G_{MD} [kcal/mol]')
ylabel('N molecules')

title('VdW/SPC/HNC')
NL=sprintf('\n');
text(min(X),max(H),['STD = ' num2str(S) ' kcal/mol' NL 'MEAN = ' num2str(M) ' kcal/mol']);

print -dpng VdW_hyst.png

Chandler = doFit('No','molecules.lst','MDFT_Hydro_Both.chain',',hydro_type=C,closure=HNC,water_model=spc','water_model,closure,hydro_type');


recalc_nocharge= doFit('No','molecules.lst','Recalc_HS.chain');
plotFit(recalc_nocharge,'RecalcNoCharge','MD')

recalc_charge= doFit('Yes','molecules.lst','Recalc_HS.chain');
plotFit(recalc_charge,'RecalcWithCharge','MD')


header ={ 'method','RMSD','a','b','N diverged (of 504 )' };

cname={'','','With Charges','',''};
C = [ recalc_charge.header celificate([recalc_charge.Rmsd  recalc_charge.a recalc_charge.b  sum(isnan(recalc_charge.FE))'])];
ncname={'','','No Charges','',''};
NC = [ recalc_nocharge.header celificate([recalc_nocharge.Rmsd  recalc_nocharge.a recalc_nocharge.b  sum(isnan(recalc_nocharge.FE))'])];

cell2tex([cname;header;C;ncname;header;NC],'RecalcNoHydro.tex',true);



chargesHS = doFit('Yes','molecules.lst','MDFT_NoHydro.chain',',hs=T');
nochargesHS = doFit('No','molecules.lst','MDFT_NoHydro.chain',',hs=T');

plotFit(nochargeHS,'HS_NoCharge','MD')
plotFit(chargesHS,'HS_WithCharge','MD')

AC=cell2mat(charges.a_list);
charges.a = AC(1:2:end);
charges.b = AC(2:2:end);


ANC = cell2mat(nocharges.a_list);
nocharges.a = ANC(1:2:end);
nocharges.b = ANC(2:2:end);


header ={ 'method','RMSD','a','b','N diverged (of 504 )' };

cname={'','','With Charges','',''};
C = [ charges.header celificate([charges.Rmsd  charges.a charges.b  sum(isnan(charges.FE))'])];
ncname={'','','No Charges','',''};
NC = [ nocharges.header celificate([nocharges.Rmsd  nocharges.a nocharges.b  sum(isnan(nocharges.FE))'])];


cell2tex([cname;header;C;ncname;header;NC],'NoHydro.tex',true);




H =  [ header celificate([hydro.Rmsd  hydro.a hydro.b  sum(isnan(hydro.FE))'])];


% ===============   Hydrophobicity ================== 



figure
plotPMV(hydro.rhoV,hydro.ExpN,hydro.FE,hydro.header,'MD VdW');
title('No Charges, No Correction')
print -dpng pmv_hydro_nocorr.png

H =  [ hydro.header celificate([hydro.Rmsd  hydro.a hydro.b  sum(isnan(hydro.FE))'])];

hname={'','','Hydrophobisity, No charge','',''};
hydroheader ={ 'method','RMSD','a','b','N diverged (of 50 )' };

cell2tex([cname;header;C;ncname;header;NC;hname;hydroheader;H],'Results.tex');



