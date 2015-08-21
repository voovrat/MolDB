[DU,EU]=unit_names;


RefC = moldb_getRef('clusters.chain');
RefS = moldb_getRef('sinja.chain');

Ref = [RefC; RefS];

E = cell_get_values(moldb_load(Ref,'energy.dat'));
Etot = cell_get_values(moldb_load(Ref,'energy_tot.dat'));

%%%%%%%% hydronium

Inw = moldb_filter(Ref,{'prog','nwchem'});
Iab = moldb_filter(Ref,{'prog','abinit'});
Ilm = moldb_filter(Ref,{'prog','lammps'});
Iru = moldb_filter(Ref,{'prog','runner'});
Ifh = moldb_filter(Ref,{'prog','fhaims'});
Ian = moldb_filter(Ref,{'prog','abinit_neutral'});


Ip = moldb_filter(Ref,{'proton','T'});
In = moldb_filter(Ref,{'proton','F'});

Inwp = intersect(Inw,Ip);
Inwn = intersect(Inw,In);

Iabp = intersect(Iab,Ip);
Iabn = intersect(Iab,In);

Ilmp = intersect(Ilm,Ip);
Ilmn = intersect(Ilm,In);

Irup = intersect(Iru,Ip);
Irun = intersect(Iru,In);

Ifhp = intersect(Ifh,Ip);
Ifhn = intersect(Ifh,In);

Ianp = intersect(Ian,Ip);
Iann = intersect(Ian,In);


dE_nw0 = load('E_nwchem_c1p.dat') - load('E_nwchem_c1n.dat');
dE_ab0 = load('E_abinit_c1p.dat') - load('E_abinit_c1n.dat');
dE_lm0 = load('E_lammps_c1p.dat') - load('E_lammps_c1n.dat');
dE_lt0 = load('E_lammpstot_c1p.dat') - load('E_lammpstot_c1n.dat');

dE_an0 = load('E_abinit_neutral_c1p.dat') - load('E_abinit_c1n.dat');

dEnn_nw = cell(7,1);
dEnn_ab = cell(7,1);
dEnn_lm = cell(7,1);
dEnn_lt = cell(7,1);
dEnn_an = cell(7,1);

for i=1:length(Inwp)
   
   cc = Ref{Inwp(i)}.parameters.molecule
   csize = length(cc)/2;

   % nwchem
   
   ii_nwp = moldb_filter(Ref(Inwp),{'molecule',cc});
   ii_nwn = moldb_filter(Ref(Inwn),{'molecule',cc});

   dE = E(Inwp(ii_nwp)) - E(Inwn(ii_nwn)) - dE_nw0;
   dEnn_nw{csize} = [ dEnn_nw{csize} dE ];


   %abinit

   ii_abp = moldb_filter(Ref(Iabp),{'molecule',cc});
   ii_abn = moldb_filter(Ref(Iabn),{'molecule',cc});
 
   dE = E(Iabp(ii_abp)) - E(Iabn(ii_abn)) - dE_ab0;
   dEnn_ab{csize} = [ dEnn_ab{csize} dE ];

   %lammps

   ii_lmp = moldb_filter(Ref(Ilmp),{'molecule',cc});
   ii_lmn = moldb_filter(Ref(Ilmn),{'molecule',cc});
 
   dE = E(Ilmp(ii_lmp)) - E(Ilmn(ii_lmn)) - dE_lm0;
   dEnn_lm{csize} = [ dEnn_lm{csize} dE ];

   %lammpstot

   ii_lmp = moldb_filter(Ref(Ilmp),{'molecule',cc});
   ii_lmn = moldb_filter(Ref(Ilmn),{'molecule',cc});
 
   dE = Etot(Ilmp(ii_lmp)) - Etot(Ilmn(ii_lmn)) - dE_lt0;
   dEnn_lt{csize} = [ dEnn_lt{csize} dE ];

   % abinit_neutral
   ii_anp = moldb_filter(Ref(Ianp),{'molecule',cc});
   ii_ann = moldb_filter(Ref(Iann),{'molecule',cc});
 
   dE = E(Ianp(ii_anp)) - E(Iann(ii_ann)) - dE_an0;
   dEnn_an{csize} = [ dEnn_an{csize} dE ];



end

dEn_nw = zeros(7,1);
dEn_ab = zeros(7,1);
dEn_lm = zeros(7,1);
dEn_lt = zeros(7,1);
dEn_an = zeros(7,1);

for i=1:7
  dEn_nw(i) = sum( dEnn_nw{i}) / Cmn(i-1,6);
  dEn_ab(i) = sum( dEnn_ab{i}) / Cmn(i-1,6);
  dEn_lm(i) = sum( dEnn_lm{i}) / Cmn(i-1,6);
  dEn_lt(i) = sum( dEnn_lt{i}) / Cmn(i-1,6);
  dEn_an(i) = sum( dEnn_an{i}) / Cmn(i-1,6);

end


dEn_nw(4) = (dEn_nw(3) + dEn_nw(5))/2;

plot([dEn_nw dEn_ab dEn_lm dEn_lt dEn_an]*unit2unit(EU,'eV','kcal/mol'))



ddE_ab = dEn_ab(end)-dEn_ab(1)
ddE_nw = dEn_nw(end) - dEn_nw(1)
ddE_lm = dEn_lm(end)-dEn_lm(1)
ddE_lt = dEn_lt(end) - dEn_lt(1)
ddE_an = dEn_an(end)-dEn_an(1)


[DU,EU]=unit_names;

ddE_ab*unit2unit(EU,'eV','kcal/mol')
ddE_nw*unit2unit(EU,'eV','kcal/mol')
ddE_lm*unit2unit(EU,'eV','kcal/mol')
ddE_lt*unit2unit(EU,'eV','kcal/mol')
ddE_an*unit2unit(EU,'eV','kcal/mol')


plot([dEn_nw-dEn_nw(end) dEn_ab-dEn_ab(end) dEn_lm-dEn_lm(end) dEn_an-dEn_an(end)]*unit2unit(EU,'eV','kcal/mol'))

%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    dimer cluster

cname=mycellfun(@(x)x.parameters.molecule,Ref)

II =  find(...
       cell_get_values(...
           mycellfun(@(x)length(x),strfind(cname,'c0'),0)...
       )...
      );

RRef = Ref(II);
EE=E(II);

IInw = moldb_filter(RRef,{'prog','nwchem'});
IIab = moldb_filter(RRef,{'prog','abinit'});
IIlm = moldb_filter(RRef,{'prog','lammps'});
IIfh = moldb_filter(RRef,{'prog','fhaims'});
IIru = moldb_filter(RRef,{'prog','runner'});
IIan = moldb_filter(RRef,{'prog','abinit_neutral'});


IIp = moldb_filter(RRef,{'proton','T'});
IIn = moldb_filter(RRef,{'proton','F'});

IInwp = intersect(IInw,IIp);
IInwn = intersect(IInw,IIn);

IIabp = intersect(IIab,IIp);
IIabn = intersect(IIab,IIn);

IIlmp = intersect(IIlm,IIp);
IIlmn = intersect(IIlm,IIn);

IIfhp = intersect(IIfh,IIp);
IIfhn = intersect(IIfh,IIn);

IIrup = intersect(IIru,IIp);
IIrun = intersect(IIru,IIn);

IIanp = intersect(IIan,IIp);
IIann = intersect(IIan,IIn);


ii0_nw_p = IInwp(moldb_filter(RRef(IInwp),{'molecule','c0'}))
ii0_nw_n = IInwn(moldb_filter(RRef(IInwn),{'molecule','c0'}))

ii0_ab_p = IIabp(moldb_filter(RRef(IIabp),{'molecule','c0'}))
ii0_ab_n = IIabn(moldb_filter(RRef(IIabn),{'molecule','c0'}))

ii0_lm_p = IIlmp(moldb_filter(RRef(IIlmp),{'molecule','c0'}))
ii0_lm_n = IIlmn(moldb_filter(RRef(IIlmn),{'molecule','c0'}))


ii0_fh_p = IIfhp(moldb_filter(RRef(IIfhp),{'molecule','c0'}))
ii0_fh_n = IIfhn(moldb_filter(RRef(IIfhn),{'molecule','c0'}))

ii0_ru_p = IIrup(moldb_filter(RRef(IIrup),{'molecule','c0'}))
ii0_ru_n = IIrun(moldb_filter(RRef(IIrun),{'molecule','c0'}))

ii0_an_p = IIanp(moldb_filter(RRef(IIanp),{'molecule','c0'}))
ii0_an_n = IIann(moldb_filter(RRef(IIann),{'molecule','c0'}))


dEE_nw0 = EE(ii0_nw_p) - EE(ii0_nw_n);
dEE_ab0 = EE(ii0_ab_p) - EE(ii0_ab_n);
dEE_lm0 = EE(ii0_lm_p) - EE(ii0_lm_n) ;
dEE_fh0 = EE(ii0_fh_p) - EE(ii0_fh_n);
dEE_ru0 = EE(ii0_ru_p) - EE(ii0_ru_n) ;
dEE_an0 = EE(ii0_an_p) - EE(ii0_an_n) ;


%dEE_lt0 = load('E_lammpstot_c1p.dat') - load('E_lammpstot_c1n.dat');



dEEnn_nw = cell(6,1);
dEEnn_ab = cell(6,1);
dEEnn_lm = cell(6,1);
dEEnn_ru = cell(6,1);
dEEnn_fh = cell(6,1);
dEEnn_an = cell(6,1);

%dEEnn_lt = cell(7,1);

for i=1:length(IInwp)
   
   cc = RRef{IInwp(i)}.parameters.molecule
   if strcmp(cc,'c0') 
      continue
   end

   csize = length(cc)/2 - 1;

   % nwchem
   
   ii_nwp = moldb_filter(RRef(IInwp),{'molecule',cc});
   ii_nwn = moldb_filter(RRef(IInwn),{'molecule',cc});

   dEE = EE(IInwp(ii_nwp)) - EE(IInwn(ii_nwn)) - dEE_nw0;
   dEEnn_nw{csize} = [ dEEnn_nw{csize} dEE ];


   %abinit

   ii_abp = moldb_filter(RRef(IIabp),{'molecule',cc});
   ii_abn = moldb_filter(RRef(IIabn),{'molecule',cc});
 
   dEE = EE(IIabp(ii_abp)) - EE(IIabn(ii_abn)) - dEE_ab0;
   dEEnn_ab{csize} = [ dEEnn_ab{csize} dEE ];

   %lammps

   ii_lmp = moldb_filter(RRef(IIlmp),{'molecule',cc});
   ii_lmn = moldb_filter(RRef(IIlmn),{'molecule',cc});
 
   dEE = EE(IIlmp(ii_lmp)) - EE(IIlmn(ii_lmn)) - dEE_lm0;
   dEEnn_lm{csize} = [ dEEnn_lm{csize} dEE ];

  %fhaims

   ii_fhp = moldb_filter(RRef(IIfhp),{'molecule',cc});
   ii_fhn = moldb_filter(RRef(IIfhn),{'molecule',cc});
 
   dEE = EE(IIfhp(ii_fhp)) - EE(IIfhn(ii_fhn)) - dEE_fh0;
   dEEnn_fh{csize} = [ dEEnn_fh{csize} dEE ];

  %runner

   ii_rup = moldb_filter(RRef(IIrup),{'molecule',cc});
   ii_run = moldb_filter(RRef(IIrun),{'molecule',cc});
 
   dEE = EE(IIrup(ii_rup)) - EE(IIrun(ii_run)) - dEE_ru0;
   dEEnn_ru{csize} = [ dEEnn_ru{csize} dEE ];


   % abinit_neutral
   ii_anp = moldb_filter(RRef(IIanp),{'molecule',cc});
   ii_ann = moldb_filter(RRef(IIann),{'molecule',cc});
 
   dEE = EE(IIanp(ii_anp)) - EE(IIann(ii_ann)) - dEE_an0;
   dEEnn_an{csize} = [ dEEnn_an{csize} dEE ];


   %lammpstot

%   ii_lmp = moldb_filter(Ref(Ilmp),{'molecule',cc});
%   ii_lmn = moldb_filter(Ref(Ilmn),{'molecule',cc});
 
%   dE = Etot(Ilmp(ii_lmp)) - Etot(Ilmn(ii_lmn)) - dE_lt0;
%   dEnn_lt{csize} = [ dEnn_lt{csize} dE ];



end




dEEn_nw = zeros(6,1);
dEEn_ab = zeros(6,1);
dEEn_lm = zeros(6,1);
dEEn_fh = zeros(6,1);
dEEn_ru = zeros(6,1);
dEEn_an = zeros(6,1);

for i=1:6
  dEEn_nw(i) = sum( dEEnn_nw{i}) / Cmn(i-1,5);
  dEEn_ab(i) = sum( dEEnn_ab{i}) / Cmn(i-1,5);
  dEEn_lm(i) = sum( dEEnn_lm{i}) / Cmn(i-1,5);
  dEEn_fh(i) = sum( dEEnn_fh{i}) / Cmn(i-1,5);
  dEEn_ru(i) = sum( dEEnn_ru{i}) / Cmn(i-1,5);
  dEEn_an(i) = sum( dEEnn_an{i}) / Cmn(i-1,5);
%  dEEn_lt(i) = sum( dEnn_lt{i}) / Cmn(i-1,6);
end


dEEn_nw(3) = (dEEn_nw(2) + dEEn_nw(4))/2;


plot([dEEn_nw dEEn_ab dEEn_lm dEEn_fh dEEn_ru dEEn_an]*unit2unit(EU,'eV','kcal/mol'))


%%%%%%%

ddEE_ab = dEEn_ab(end) - dEEn_ab(1)
ddEE_nw = dEEn_nw(end) - dEEn_nw(1)
ddEE_lm = dEEn_lm(end) - dEEn_lm(1)
ddEE_ru = dEEn_ru(end) - dEEn_ru(1)
ddEE_fh = dEEn_fh(end) - dEEn_fh(1)
ddEE_an = dEEn_an(end) - dEEn_an(1);
%ddE_lt = dEn_lt(end) - dEn_lt(1)


[DU,EU]=unit_names;

ddEE_ab*unit2unit(EU,'eV','kcal/mol')
ddEE_nw*unit2unit(EU,'eV','kcal/mol')
ddEE_lm*unit2unit(EU,'eV','kcal/mol')
ddEE_ru*unit2unit(EU,'eV','kcal/mol')
ddEE_fh*unit2unit(EU,'eV','kcal/mol')
ddEE_an*unit2unit(EU,'eV','kcal/mol')
%ddE_lt*unit2unit(EU,'eV','kcal/mol')


%D1 = [dEn_nw-dEn_nw(end) dEn_ab-dEn_ab(end) dEn_lm-dEn_lm(end) ] *unit2unit(EU,'eV','kcal/mol');
%D2 = [dEEn_nw-dEEn_nw(end) dEEn_ab-dEEn_ab(end) dEEn_lm-dEEn_lm(end)  dEEn_fh-dEEn_fh(end) ] * unit2unit(EU,'eV','kcal/mol');

%plot((1:7),D1,(2:7),D2)


%plot([dEn_nw-dEn_nw(end) dEn_ab-dEn_ab(end) dEn_lm-dEn_lm(end)]*unit2unit(EU,'eV','kcal/mol'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%% trimer cluster


cname=mycellfun(@(x)x.parameters.molecule,Ref)

III =  find(...
       cell_get_values(...
           mycellfun(@(x)length(x),strfind(cname,'c0c1'),0)...
       )...
      );

RRRef = Ref(III);
EEE=E(III);

IIInw = moldb_filter(RRRef,{'prog','nwchem'});
IIIab = moldb_filter(RRRef,{'prog','abinit'});
IIIlm = moldb_filter(RRRef,{'prog','lammps'});
IIIfh = moldb_filter(RRRef,{'prog','fhaims'});
IIIru = moldb_filter(RRRef,{'prog','runner'});
IIIan = moldb_filter(RRRef,{'prog','abinit_neutral'});

IIIp = moldb_filter(RRRef,{'proton','T'});
IIIn = moldb_filter(RRRef,{'proton','F'});

IIInwp = intersect(IIInw,IIIp);
IIInwn = intersect(IIInw,IIIn);

IIIabp = intersect(IIIab,IIIp);
IIIabn = intersect(IIIab,IIIn);

IIIlmp = intersect(IIIlm,IIIp);
IIIlmn = intersect(IIIlm,IIIn);

IIIfhp = intersect(IIIfh,IIIp);
IIIfhn = intersect(IIIfh,IIIn);

IIIrup = intersect(IIIru,IIIp);
IIIrun = intersect(IIIru,IIIn);

IIIanp = intersect(IIIan,IIIp);
IIIann = intersect(IIIan,IIIn);


iii0_nw_p = IIInwp(moldb_filter(RRRef(IIInwp),{'molecule','c0c1'}))
iii0_nw_n = IIInwn(moldb_filter(RRRef(IIInwn),{'molecule','c0c1'}))

iii0_ab_p = IIIabp(moldb_filter(RRRef(IIIabp),{'molecule','c0c1'}))
iii0_ab_n = IIIabn(moldb_filter(RRRef(IIIabn),{'molecule','c0c1'}))

iii0_lm_p = IIIlmp(moldb_filter(RRRef(IIIlmp),{'molecule','c0c1'}))
iii0_lm_n = IIIlmn(moldb_filter(RRRef(IIIlmn),{'molecule','c0c1'}))

iii0_fh_p = IIIfhp(moldb_filter(RRRef(IIIfhp),{'molecule','c0c1'}))
iii0_fh_n = IIIfhn(moldb_filter(RRRef(IIIfhn),{'molecule','c0c1'}))

iii0_ru_p = IIIrup(moldb_filter(RRRef(IIIrup),{'molecule','c0c1'}))
iii0_ru_n = IIIrun(moldb_filter(RRRef(IIIrun),{'molecule','c0c1'}))

iii0_an_p = IIIanp(moldb_filter(RRRef(IIIanp),{'molecule','c0c1'}))
iii0_an_n = IIIann(moldb_filter(RRRef(IIIann),{'molecule','c0c1'}))


dEEE_nw0 = EEE(iii0_nw_p) - EEE(iii0_nw_n);
dEEE_ab0 = EEE(iii0_ab_p) - EEE(iii0_ab_n);
dEEE_lm0 = EEE(iii0_lm_p) - EEE(iii0_lm_n) ;
dEEE_fh0 = EEE(iii0_fh_p) - EEE(iii0_fh_n);
dEEE_ru0 = EEE(iii0_ru_p) - EEE(iii0_ru_n) ;
dEEE_an0 = EEE(iii0_an_p) - EEE(iii0_an_n) ;

%dEE_lt0 = load('E_lammpstot_c1p.dat') - load('E_lammpstot_c1n.dat');



dEEEnn_nw = cell(5,1);
dEEEnn_ab = cell(5,1);
dEEEnn_lm = cell(5,1);
dEEEnn_ru = cell(5,1);
dEEEnn_fh = cell(5,1);
dEEEnn_an = cell(5,1);


%dEEnn_lt = cell(7,1);

for i=1:length(IIInwp)
   
   cc = RRRef{IIInwp(i)}.parameters.molecule
   if strcmp(cc,'c0c1') 
      continue
   end

   csize = length(cc)/2 - 2;

   % nwchem
   
   ii_nwp = moldb_filter(RRRef(IIInwp),{'molecule',cc});
   ii_nwn = moldb_filter(RRRef(IIInwn),{'molecule',cc});

   dEEE = EEE(IIInwp(ii_nwp)) - EEE(IIInwn(ii_nwn)) - dEEE_nw0;
   dEEEnn_nw{csize} = [ dEEEnn_nw{csize} dEEE ];


   %abinit

   ii_abp = moldb_filter(RRRef(IIIabp),{'molecule',cc});
   ii_abn = moldb_filter(RRRef(IIIabn),{'molecule',cc});
 
   dEEE = EEE(IIIabp(ii_abp)) - EEE(IIIabn(ii_abn)) - dEEE_ab0;
   dEEEnn_ab{csize} = [ dEEEnn_ab{csize} dEEE ];

   %lammps

   ii_lmp = moldb_filter(RRRef(IIIlmp),{'molecule',cc});
   ii_lmn = moldb_filter(RRRef(IIIlmn),{'molecule',cc});
 
   dEEE = EEE(IIIlmp(ii_lmp)) - EEE(IIIlmn(ii_lmn)) - dEEE_lm0;
   dEEEnn_lm{csize} = [ dEEEnn_lm{csize} dEEE ];

   %fhaims

   ii_fhp = moldb_filter(RRRef(IIIfhp),{'molecule',cc});
   ii_fhn = moldb_filter(RRRef(IIIfhn),{'molecule',cc});
 
   dEEE = EEE(IIIfhp(ii_fhp)) - EEE(IIIfhn(ii_fhn)) - dEEE_fh0;
   dEEEnn_fh{csize} = [ dEEEnn_fh{csize} dEEE ];

   %runner

   ii_rup = moldb_filter(RRRef(IIIrup),{'molecule',cc});
   ii_run = moldb_filter(RRRef(IIIrun),{'molecule',cc});
 
   dEEE = EEE(IIIrup(ii_rup)) - EEE(IIIrun(ii_run)) - dEEE_ru0;
   dEEEnn_ru{csize} = [ dEEEnn_ru{csize} dEEE ];


   %abinit_neutral

   ii_anp = moldb_filter(RRRef(IIIanp),{'molecule',cc});
   ii_ann = moldb_filter(RRRef(IIIann),{'molecule',cc});
 
   dEEE = EEE(IIIanp(ii_anp)) - EEE(IIIann(ii_ann)) - dEEE_an0;
   dEEEnn_an{csize} = [ dEEEnn_an{csize} dEEE ];


   %lammpstot



%   ii_lmp = moldb_filter(Ref(Ilmp),{'molecule',cc});
%   ii_lmn = moldb_filter(Ref(Ilmn),{'molecule',cc});
 
%   dE = Etot(Ilmp(ii_lmp)) - Etot(Ilmn(ii_lmn)) - dE_lt0;
%   dEnn_lt{csize} = [ dEnn_lt{csize} dE ];



end




dEEEn_nw = zeros(5,1);
dEEEn_ab = zeros(5,1);
dEEEn_lm = zeros(5,1);
dEEEn_fh = zeros(5,1);
dEEEn_ru = zeros(5,1);
dEEEn_an = zeros(5,1);

%dEEn_lt = zeros(7,1);

for i=1:5
  dEEEn_nw(i) = sum( dEEEnn_nw{i}) / Cmn(i-1,4);
  dEEEn_ab(i) = sum( dEEEnn_ab{i}) / Cmn(i-1,4);
  dEEEn_lm(i) = sum( dEEEnn_lm{i}) / Cmn(i-1,4);
  dEEEn_fh(i) = sum( dEEEnn_fh{i}) / Cmn(i-1,4);
  dEEEn_ru(i) = sum( dEEEnn_ru{i}) / Cmn(i-1,4);
  dEEEn_an(i) = sum( dEEEnn_an{i}) / Cmn(i-1,4);

%  dEEn_lt(i) = sum( dEnn_lt{i}) / Cmn(i-1,6);
end


dEEEn_nw(2) = (dEEEn_nw(1) + dEEEn_nw(3))/2;


plot([dEEEn_nw dEEEn_ab dEEEn_lm dEEEn_fh dEEEn_an]*unit2unit(EU,'eV','kcal/mol'))


%%%%%%%

ddEEE_ab = dEEEn_ab(end) - dEEEn_ab(1)
ddEEE_nw = dEEEn_nw(end) - dEEEn_nw(1)
ddEEE_lm = dEEEn_lm(end) - dEEEn_lm(1)
ddEEE_ru = dEEEn_ru(end) - dEEEn_ru(1)
ddEEE_fh = dEEEn_fh(end) - dEEEn_fh(1)
ddEEE_an = dEEEn_an(end) - dEEEn_an(1)

%ddE_lt = dEn_lt(end) - dEn_lt(1)


[DU,EU]=unit_names;

ddEEE_ab*unit2unit(EU,'eV','kcal/mol')
ddEEE_nw*unit2unit(EU,'eV','kcal/mol')
ddEEE_lm*unit2unit(EU,'eV','kcal/mol')
ddEEE_ru*unit2unit(EU,'eV','kcal/mol')
ddEEE_fh*unit2unit(EU,'eV','kcal/mol')
ddEEE_an*unit2unit(EU,'eV','kcal/mol')


%ddE_lt*unit2unit(EU,'eV','kcal/mol')


D1 = [dEn_nw-dEn_nw(end) dEn_ab-dEn_ab(end) dEn_lm-dEn_lm(end)] *unit2unit(EU,'eV','kcal/mol');
D2 = [dEEn_nw-dEEn_nw(end) dEEn_ab-dEEn_ab(end) dEEn_lm-dEEn_lm(end) dEEn_fh-dEEn_fh(end)] * unit2unit(EU,'eV','kcal/mol');

D3 = [dEEEn_nw-dEEEn_nw(end) dEEEn_ab-dEEEn_ab(end) dEEEn_lm-dEEEn_lm(end) dEEE_fh-dEEE_fh(end)] * unit2unit(EU,'eV','kcal/mol');


plot((1:7),D1,(2:7),D2,(3:7),D3)


KK = unit2unit(EU,'eV','kcal/mol');

%%%%%%%%%%%%%%%%

head = {'m','nwchem','abinit','lammps'}
A1 = [ dEn_nw dEn_ab dEn_lm dEn_an ] * KK;
X1 = (1:7)';
O1 = ones(size(X1));
B1 = A1 - O1*A1(end,:);

%C1 = [ head; celificate(A1) ];
C1 = [ celificate([X1 A1 ]) ; {'max-min'} celificate(max(A1) - min(A1)) ];

print_cell(C1,'c1p.tex','&','\\')


%%%%%%%%%%%%%%%%

head = {'m','nwchem','abinit','lammps','fhaims','abinit_neutral'}
A2 = [ dEEn_nw dEEn_ab dEEn_lm dEEn_fh dEEn_an ] * KK;
X2 = (1:6)';
O2 = ones(size(X2));
B2 = A2 - O2*A2(end,:);

%C1 = [ head; celificate(A1) ];
C2 = celificate([X2 A2 ]);
C2 = [ celificate([X2 A2 ]) ; {'max-min'} celificate(max(A2) - min(A2)) ];


print_cell(C2,'c2p.tex','&','\\')

%%%%%%%%%%%%%%%%

head = {'m','nwchem','abinit','lammps'}
A3 = [ dEEEn_nw dEEEn_ab dEEEn_lm dEEEn_fh dEEEn_an ] * KK;
X3 = (1:5)';
O3 = ones(size(X3));
B3 = A3 - O3*A3(end,:);

%C1 = [ head; celificate(A1) ];
C3 = celificate([X3 A3]);
C3 = [ celificate([X3 A3 ]) ; {'max-min'} celificate(max(A3) - min(A3)) ];


print_cell(C3,'c3p.tex','&','\\')


%%%%%%%%%

kT = unit2unit(EU,'K','kcal/mol')*298.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hnd = plot((1:7),B1,[1 7],[kT kT],[1 7],[-kT -kT])
xlabel('m (bulk cluster size)')
ylabel('$\Delta E_{meth}^m - \Delta E_{meth}^7$ [kcal/mol]')
legend('nwchem','abinit','lammps(NN)','abinit-neutral','kT')

set(hnd(1:4),'LineWidth',3)

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','^')
set(hnd(4),'Marker','x','Color',[1 0 1])
set(hnd([5 6]),'LineStyle','--','Color',[0 0 0])

printpng('fig_c1p.png')

%print -F:16 c1p.png
%print -F:16 -depslatex  c1p.eps
%print -F:16 -depslatexstandalone  c1p.tex
%system('latex c1p.tex')
%system('dvipdf c1p.dvi')
%system('evince c1p.pdf')
%system('convert -density 150 c1p.pdf c1p_X.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hnd = plot((1:6),B2,[1 6],[kT kT],[1 6],[-kT -kT])
xlabel('m (bulk cluster size)')
ylabel('$\Delta E_{meth}^m - \Delta E_{meth}^6$ [kcal/mol]')
legend('nwchem','abinit','lammps(RuNNer)','FHaims','abinit-neutral','kT')

set(hnd(1:5),'LineWidth',3)

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','^')
set(hnd(4),'Marker','d')
set(hnd(5),'Marker','x','Color',[1 0 1])
set(hnd([ 6 7]),'LineStyle','--','Color',[0 0 0])
legend('boxoff')

printpng('fig_c2p.png');

%print -F:16 c2p.png

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hnd = plot((1:5),B3,[1 5],[kT kT],[1 5],[-kT -kT])
xlabel('m (bulk cluster size)')
ylabel('$\Delta E_{meth}^m - \Delta E_{meth}^5$ [kcal/mol]')
legend('nwchem','abinit','lammps(NN)','FHaims','abinit-neutral','kT')

set(hnd(1:5),'LineWidth',3)

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','^')
set(hnd(4),'Marker','d')
set(hnd(5),'Marker','x','Color',[1 0 1])
set(hnd([ 6 7]),'LineStyle','--','Color',[0 0 0])
legend('boxoff')
ylim([-6 10])
printpng('fig_c3p.png');

%print -F:16 c3p.png

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hnd = plot(B2)
xlabel('m (bulk cluster size)')
ylabel('\Delta E_{meth}^m - \Delta E_{meth}^6 [kcal/mol]')
legend('nwchem','abinit','lammps(NN)')

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','^')
set(hnd,'LineWidth',3)

print -F:16 c2p.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hnd = plot(B3)
xlabel('m (bulk cluster size)')
ylabel('\Delta E_{meth}^m - \Delta E_{meth}^5 [kcal/mol]')
legend('nwchem','abinit','lammps(NN)')

set(hnd(1),'Marker','s')
set(hnd(2),'Marker','o')
set(hnd(3),'Marker','^')
set(hnd,'LineWidth',3)

print -F:16 c3p.png





