Ref0 = moldb_getRef('clusters.chain');
Ref1 = moldb_getRef('sinja.chain');

Ref = [Ref0;Ref1];

E = cell_get_values(moldb_load(Ref,'energy.dat'));
MM = cell_get_values(moldb_load(Ref,'max_min_dist.dat'));

Inw = moldb_filter(Ref,{'prog','nwchem'});
Iab = moldb_filter(Ref,{'prog','abinit'});
Iru = moldb_filter(Ref,{'prog','runner'});
Ifh = moldb_filter(Ref,{'prog','fhaims'});
Ian = moldb_filter(Ref,{'prog','abinit_neutral'});

Ip = moldb_filter(Ref,{'proton','T'});
In = moldb_filter(Ref,{'proton','F'});

Inwp = intersect(Inw,Ip);
Iabp = intersect(Iab,Ip);
Irup = intersect(Iru,Ip);
Ifhp = intersect(Ifh,Ip);
Ianp = intersect(Ian,Ip);

Inwn = intersect(Inw,In);
Iabn = intersect(Iab,In);
Irun = intersect(Iru,In);
Ifhn = intersect(Ifh,In);
Iann = intersect(Ian,In);

Enw0 = load('E_nwchem_c1n.dat');
Eab0 = load('E_abinit_c1n.dat');
Eru0 = load('E_lammps_c1n.dat');
Efh0 = Eru0;
Ean0 = Eab0;

Nmol = cell_get_values(mycellfun(@(x)length(x.parameters.molecule)/2,Ref)) + 1;

E0 = zeros(size(E));
E0(Inw) = Enw0;
E0(Iab) = Eab0;
E0(Iru) = Eru0;
E0(Ifh) = Efh0;
E0(Ian) = Ean0;

I = [Inw; Iab; Iru; Ifh; Ian];

E1 = E - Nmol.*E0;

Ifive = find(Nmol==5);

E1( intersect(Inw,Ifive) ) = nan;


hnd=plot([E1(Inw) E1(Iab) E1(Iru) E1(Ifh) E1(Ian) ],'.')

xlabel('\textbf{cluster}');
ylabel('\textbf{normalized energy [eV]}')

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','x')
set(hnd(4),'Marker','+')
set(hnd(5),'Marker','^')

legend('nwchem','abinit','runner','FHaims','abinit-neutral','location','eastoutside')

xlim([0 256])

printpng('normalized_energy.png');

%% *********************** %%

dEnw = E1(Inwp) - E1(Inwn);
dEab = E1(Iabp) - E1(Iabn);
dEru = E1(Irup) - E1(Irun);
dEfh = E1(Ifhp) - E1(Ifhn);
dEan = E1(Ianp) - E1(Iann);

hnd=plot([dEnw dEab dEru dEfh dEan],'.')


xlabel('\textbf{cluster}');
ylabel('$E^+ - E^0$ [eV]')

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','x')
set(hnd(4),'Marker','+')
set(hnd(5),'Marker','^')

legend('nwchem','abinit','runner','FHaims','abinit-neutral','location','eastoutside')

xlim([0 128])

printpng('dE.png');

%% *********************** %%

molname = mycellfun(@(x)x.parameters.molecule,Ref(Inwp));
ii = find(strcmp(molname,'c0'));

II0 = find( cell_get_values( mycellfun (@(x)length(strfind(x,'c0')),molname))  );

II=setdiff(II0,ii);



dEEnw = dEnw(II) - dEnw(ii);
dEEab = dEab(II) - dEab(ii);
dEEru = dEru(II) - dEru(ii);
dEEfh = dEfh(II) - dEfh(ii);
dEEan = dEan(II) - dEan(ii);


hnd=plot([dEEnw dEEab dEEru dEEfh dEEan],'.')

xlabel('\textbf{cluster}');
ylabel('$\Delta E_{clust} = E^+_{clust} - E^0_{clust} - E^+_{H_5O_2+} + E^0_{H_4O_2}$ [eV]')

set(hnd(1),'Marker','s') % nw
set(hnd(2),'Marker','o') % ab
set(hnd(3),'Marker','x') % ru 
set(hnd(4),'Marker','+') % fh
set(hnd(5),'Marker','^') % an

legend('nwchem','abinit','runner','FHaims','abinit-neutral','location','eastoutside')

xlim([0 64])

printpng('dEE.png');


hnd=plot(dEEfh, [dEEfh dEEru  dEEan  ],'.')

xlabel('$\Delta E_{clust}^{FHaims}$ [eV]');
ylabel('$\Delta E_{clust}^{meth}$ [eV]')

set(hnd(1),'Marker','+') % fh
set(hnd(2),'Marker','x') % ru 
set(hnd(3),'Marker','^') % an

legend('FHaims','runner','abinit-neutral','location','eastoutside')

printpng('fh_ru_an.png')


hnd=plot(dEEfh, [dEEfh - dEEfh dEEru - dEEfh dEEan - dEEfh  ],'.')


xlabel('$\Delta E_{clust}^{FHaims}$ [eV]');
ylabel('$\Delta E_{clust}^{meth}$ [eV]')

set(hnd(1),'Marker','+') % fh
set(hnd(2),'Marker','x') % ru 
set(hnd(3),'Marker','^') % an

legend('FHaims','runner','abinit-neutral','location','eastoutside')
tinyname=mycellfun(@(x)[ '$ {}_{{}_{' x '}}$' ],molname(II));
molLabel(dEEfh,dEEru-dEEfh,tinyname)

printpng('fh_ru_an_diff.png')

%%%%%%%%%%%



mmlab=mycellfun(@(x) num2str(x),celificate(MM(Ianp(II))));
molLabel(dEEfh,dEEru-dEEfh,mmlab)


mm = MM(Ianp(II));
nnm = Nmol(Ianp(II));
%plot(mm,(dEEru - dEEf)./nnm,'.')
hnd=plot(mm,[dEEfh - dEEfh dEEru - dEEfh dEEan - dEEfh],'s')
xlabel('max(min(dist)) [\AA]')
ylabel('$\Delta E^{meth} - \Delta E^{FHaims}$ [eV]')

set(hnd(1),'Marker','+') % fh
set(hnd(2),'Marker','x') % ru 
set(hnd(3),'Marker','^') % an

legend('FHaims','runner','abinit-neutral','location','eastoutside')
%tinyname=mycellfun(@(x)[ '$ {}_{{}_{' x '}}$' ],molname(II));
%molLabel(dEEfh,dEEru-dEEfh,tinyname)

printpng('maxmindist.png')

%%%%%%%%%%%%%%%


OxD = moldb_load(Ref0(Ianp(II)),'ox_dists.dat')
DistToFirst = cell_get_values(mycellfun(@(A) min(A(3:end,1)),OxD));

hnd=plot(DistToFirst,[dEEfh - dEEfh dEEru - dEEfh dEEan - dEEfh],'s')
xlabel('Hydronium-Water distance [\AA]')
ylabel('$\Delta E^{meth} - \Delta E^{FHaims}$ [eV]')

set(hnd(1),'Marker','+') % fh
set(hnd(2),'Marker','x') % ru 
set(hnd(3),'Marker','^') % an

legend('FHaims','runner','abinit-neutral','location','eastoutside')
%tinyname=mycellfun(@(x)[ '$ {}_{{}_{' x '}}$' ],molname(II));
%molLabel(dEEfh,dEEru-dEEfh,tinyname)

printpng('neardist.png')




plot(MM(Ianp(II)),[dEEfh - dEEru] ./ Nmol(Ianp(II)),'.')

OxD = moldb_load(Ref0(Ianp(II)),'ox_dists.dat')

RedOxD = mycellfun( @(A)[0 min(A([1 2],3:end)); min(A([1 2],3:end))' A(3:end,3:end) ],OxD);

% remove zeors
RedOxD1 = mycellfun( @(A) myeval(A,'param(param==0)=inf;result=param'),RedOxD);

MaxMin = cell_get_values( mycellfun(@(A) max(min(A)),RedOxD1));




