Ref1 = moldb_getRef('clusters_new.chain');
Ref2 = moldb_getRef('clusters.chain');

I0=moldb_filter(Ref2,{'prog',{'abinit','nwchem'} });

Ref = [Ref1;Ref2(I0)];

E = cell_get_values(moldb_load(Ref,'energy.dat'));
%MM = cell_get_values(moldb_load(Ref,'max_min_dist.dat'));

Igs = moldb_filter(Ref,{'prog','gaussian'});
Iab = moldb_filter(Ref,{'prog','abinit'});
Inw = moldb_filter(Ref,{'prog','nwchem'});
Il1 = moldb_filter(Ref,{'prog','lammps_new1'});
Il2 = moldb_filter(Ref,{'prog','lammps_new2'});

Ip = moldb_filter(Ref,{'proton','T'});
In = moldb_filter(Ref,{'proton','F'});

Igsp = intersect(Igs,Ip);
Iabp = intersect(Iab,Ip);
Inwp = intersect(Inw,Ip);
Il1p = intersect(Il1,Ip);
Il2p = intersect(Il2,Ip);

Igsn = intersect(Igs,In);
Iabn = intersect(Iab,In);
Inwn = intersect(Inw,In);
Il1n = intersect(Il1,In);
Il2n = intersect(Il2,In);

ii_gs = moldb_filter(Ref,{'prog','gaussian','proton','F','molecule','c6'});

ii_ab = moldb_filter(Ref,{'prog','abinit','proton','F','molecule','c6'});

ii_nw = moldb_filter(Ref,{'prog','nwchem','proton','F','molecule','c6'});


ii_l1 = moldb_filter(Ref,{'prog','lammps_new1','proton','F','molecule','empty'});

ii_l2 = moldb_filter(Ref,{'prog','lammps_new2','proton','F','molecule','empty'});


Egs0 = E(ii_gs) / 2;
Eab0 = E(ii_ab) / 2;
Enw0 = E(ii_nw) / 2;
El10 = E(ii_l1);
El20 = E(ii_l2);


fn = @(x) iif( strcmp(x.parameters.molecule,'empty'),...
               0, ...
               length(x.parameters.molecule)/2  );


Nmol = cell_get_values(mycellfun(@(x)fn(x),Ref)) + 1;



E0 = zeros(size(E));
E0(Igs) = Egs0;
E0(Iab) = Eab0;
E0(Inw) = Enw0;
E0(Il1) = El10;
E0(Il2) = El20;


I = [Igs; Iab; Inw; Il1; Il2];

E1 = E - Nmol.*E0;

Ifive = find(Nmol==5);

E1( intersect(Inw,Ifive) ) = nan;


hnd=plot([E1(Igs) E1(Iab) E1(Inw) E1(Il1) E1(Il2) ],'.')

xlabel('\textbf{cluster}');
ylabel('\textbf{normalized energy [eV]}')

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','x')
set(hnd(4),'Marker','+')
set(hnd(5),'Marker','^')

legend('Gaussian','abinit','nwchem','NN (I)','NN(II)','location','eastoutside')

xlim([0 256])

printpng('normalized_energy_new.png');

%% *********************** %%

dEgs = E1(Igsp) - E1(Igsn);
dEab = E1(Iabp) - E1(Iabn);
dEnw = E1(Inwp) - E1(Inwn);
dEl1 = E1(Il1p) - E1(Il1n);
dEl2 = E1(Il2p) - E1(Il2n);

hnd=plot([dEgs dEab dEnw dEl1 dEl2],'.')


xlabel('\textbf{cluster}');
ylabel('$E^+ - E^0$ [eV]')

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','x')
set(hnd(4),'Marker','+')
set(hnd(5),'Marker','^')

legend('gaussian','abinit','nwchem','NN (I)','NN (II)','location','eastoutside')

xlim([0 128])

printpng('dE_new.png');

%% *********************** %%

molname = mycellfun(@(x)x.parameters.molecule,Ref(Inwp));
ii = find(strcmp(molname,'c0'));

II0 = find( cell_get_values( mycellfun (@(x)length(strfind(x,'c0')),molname))  );

II=setdiff(II0,ii);



dEEgs = dEgs(II) - dEgs(ii);
dEEab = dEab(II) - dEab(ii);
dEEnw = dEnw(II) - dEnw(ii);
dEEl1 = dEl1(II) - dEl1(ii);
dEEl2 = dEl2(II) - dEl2(ii);


hnd=plot([dEEgs dEEab dEEnw dEEl1 dEEl2],'.')

xlabel('\textbf{cluster}');
ylabel('$\Delta E_{clust} = E^+_{clust} - E^0_{clust} - E^+_{H_5O_2+} + E^0_{H_4O_2}$ [eV]')

set(hnd(1),'Marker','s') % nw
set(hnd(2),'Marker','o') % ab
set(hnd(3),'Marker','x') % ru 
set(hnd(4),'Marker','+') % fh
set(hnd(5),'Marker','^') % an

legend('gaussian','abinit','nwchem','Lammps (I)','Lammps (II)','location','eastoutside')

xlim([0 64])

printpng('dEE_new.png');


hnd=plot(dEEgs, [dEEgs dEEab dEEnw dEEl1  dEEl2  ],'.')

xlabel('$\Delta E_{clust}^{Gauss}$ [eV]');
ylabel('$\Delta E_{clust}^{meth}$ [eV]')

set(hnd(1),'Marker','.') % gs
set(hnd(1),'Marker','+') % ab
set(hnd(2),'Marker','x') % l1 
set(hnd(3),'Marker','^') % l2

legend('gaussian','abinit','nwchem','NN (I)','NN (II)','location','eastoutside')

printpng('lammps_qm_correlation_new.png')


%hnd=plot(dEEfh, [dEEfh - dEEfh dEEru - dEEfh dEEan - dEEfh  ],'.')


%xlabel('$\Delta E_{clust}^{Gauss}$ [eV]');
%ylabel('$\Delta E_{clust}^{meth}$ [eV]')

%set(hnd(1),'Marker','+') % fh
%set(hnd(2),'Marker','x') % ru 
%set(hnd(3),'Marker','^') % an

%legend('FHaims','runner','abinit-neutral','location','eastoutside')
%tinyname=mycellfun(@(x)[ '$ {}_{{}_{' x '}}$' ],molname(II));
%molLabel(dEEfh,dEEru-dEEfh,tinyname)

%printpng('lammps_QM_correlation_new.png')

