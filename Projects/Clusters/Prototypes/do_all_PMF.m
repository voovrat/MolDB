Ref = moldb_getRef('PMF1.chain');

E = cell_get_values(moldb_load(Ref,'energy.dat'));

d = cell_get_values(mycellfun(@(x)str2num(x.parameters.d),Ref));

I1 = moldb_filter(Ref,{'dist','2.5'});
I2 = moldb_filter(Ref,{'dist','2.8'});
I3 = moldb_filter(Ref,{'dist','3'});
I4 = moldb_filter(Ref,{'dist','3.2'});
I5 = moldb_filter(Ref,{'dist','3.5'});

dd = d(I1);
hnd=plot(dd,[E(I1) E(I2) E(I3) E(I4) E(I5)]);
set(hnd,'LineWidth',3);
legend('$d_{OO} =2.5$ \AA','$d_{OO} = 2.8$ A','$d_{OO} = 3 $\AA','$d_{OO}=3.2$ \AA','$d_{OO} = 3.5$ \AA');

xlabel('Proton-oxygen distancs [\AA]')
ylabel('Energy [eV]');

legend('boxoff');
legend('location','eastoutside');
printpng('proton_transfer_energy.png')

minE = min(E);


DAT = [dd E(I1) E(I2) E(I3) E(I4) E(I5) ]
save -ascii proton_transfer_energy.dat DAT

[DU,EU]=unit_names;

%ylim([minE-0.1 minE+0.1]);


%%%%%%%%%%%%%%%%

[XYZ,at] = read_xyz('oktamer_proton.xyz');
XYZn = XYZ(1:24,:);

OX = XYZn(1:3:end,:);

DOX = mydist_fast(OX)

DOX(DOX<1e-8) = inf;

DOX1 = DOX;
DOX1(DOX1==inf) = 0;
hist(DOX1(:),30);


DOX * 1.1;

SYS = zeros(0,3);


for i=1:8
  if i==8
    nhydr = 3;
  else
    nhydr = 2;
  end
 
  ox = 3*(i-1)+1; 
  MOL0 = XYZ(ox:ox+nhydr,:);
  size(MOL0)
 
  c0 = MOL0(1,:);
  c1 = c0*1.3;
  MOL = MOL0 + ones(nhydr+1,1)*(c1 - c0);
  
  SYS = [SYS; MOL];
end

write_xyz(at,SYS,'oktamer_sparse.xyz');


% 10x10x10 --> 10x10x10x0.73  sp V 1  =  33.33   V X  
%//730 * 1 = 33.33*X
%//X = 730/33.33;
%//X^(1/3) = 2.7979




