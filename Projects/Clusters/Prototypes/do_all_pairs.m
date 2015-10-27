Ref = moldb_getRef('pair.chain')

E = cell_get_values(moldb_load(Ref,'energy.dat'));

Sing = load('OX.idx');

NSing = size(Sing,1);

SingNam = cell(NSing,1);

for i = 1:NSing
   SingNam{i} = sprintf('w%0.0fw%0.0fw%0.0f',Sing(i,1),Sing(i,2),Sing(i,3));

end

N = size(Sing,1);


lmp1 = struct;
lmp1.name = 'lammps_new1';

lmp2 = struct;
lmp2.name = 'lammps_new2';

meth = {lmp1, lmp2};

for meth_id = 1:length(meth)

  meth{meth_id}.H = zeros(NSing,NSing);

  for i = 1:NSing

  
     i1p = moldb_filter(Ref,{'molecule',SingNam{i},'prog',meth{meth_id}.name,'proton','T'} );
     E1p = E(i1p);

     i1n = moldb_filter(Ref,{'molecule',SingNam{i},'prog',meth{meth_id}.name,'proton','F'} );
     E1n = E(i1n);

     dE1 = E1p - E1n;


     meth{meth_id}.H(i,i) = dE1;

     for j=i+1:NSing

     
   
        i2p = moldb_filter(Ref,{'molecule',SingNam{j},'prog',meth{meth_id}.name,'proton','T'} );
        E2p = E(i2p);

        i2n = moldb_filter(Ref,{'molecule',SingNam{j},'prog',meth{meth_id}.name,'proton','F'} );
        E2n = E(i2n);

        dE2 = E2p - E2n;


   %--------------%

        Clust1 = Sing(i,:);
        Clust2 = Sing(j,:);      
 
        Clust12 = unique( [Clust1, Clust2 ]);
        Nam12 = '';
        for ii=1:length(Clust12)
           Nam12 = [Nam12 'w' num2str(Clust12(ii)) ];
        end

   %------------%
        i12p = moldb_filter(Ref,{'molecule',Nam12,'prog',meth{meth_id}.name,'proton','T'} );
        E12p = E(i12p);

        i12n = moldb_filter(Ref,{'molecule',Nam12,'prog',meth{meth_id}.name,'proton','F'} );
        E12n = E(i12n);

        dE12 = E12p - E12n;

  %--------------%
       
       % if dE12 < min(dE1,dE2)
           H12 = sqrt( (dE12 - dE1)*(dE12-dE2) ); 
       % else
       %    H12 = 0;
       % end
        meth{meth_id}.H(i,j) = H12;
        meth{meth_id}.H(j,i) = H12;


      end
  end



end




