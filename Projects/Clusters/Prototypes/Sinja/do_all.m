

% read molecules
[A,atoms] = read_xyz('h5o2-1.xyz');
B = read_xyz('h5o2-2.xyz');

% centrate ( set 1st oxygen to the origin)
A1 = A - ones(7,1)*A(1,:)
B1 = B - ones(7,1)*B(1,:)

%aa=[atoms;atoms];
%write_xyz(aa,[A1;B1],'AB.xyz')

phi1 = atan2(A1(2,2),A1(2,1))
phi2 = atan2(B1(2,2),B1(2,1))

% rotate them both to the Ox axis

% rotate A1 by -phi1
RotPhi1 = [ cos(phi1) sin(phi1) 0 ; 
            -sin(phi1) cos(phi1)  0;
               0        0        1];

A2 = (RotPhi1 * A1' )';

% rotate B1 by -phi2
RotPhi2 = [ cos(phi2) sin(phi2) 0 ; 
            -sin(phi2) cos(phi2)  0;
               0        0        1];

B2 = (RotPhi2 * B1' )';


%write_xyz(aa,[A2;B2],'AB2.xyz');

%now y of the second oxygen of both molecules is 0
% we calculate thetas
th1 = atan2(A2(2,3),A2(2,1));
th2 = atan2(B2(2,3),B2(2,1));

% th2 --> th1
% th2 + (th1 - th2)
% rotate B2 by (th1 - th2)
RotTh= [ cos(th1-th2)  0   -sin(th1-th2);
          0            1        0
         sin(th1-th2)  0    cos(th1-th2) ];

B3 = (RotTh*B2')';

%write_xyz(aa,[A2;B3],'AB3.xyz')

% now oxygens of B3 lie on the same line as oxygens of A2
% we can rotate back B3 by phi1 to make it similar to A1


RotPh = [ cos(phi1) -sin(phi1) 0; 
          sin(phi1) cos(phi1)  0;
           0          0        1];

B4 = (RotPh*B3')';

%write_xyz(aa,[A1;B4],'AB4.xyz');

% and now we can shift back the B4 to make it coinsiding with A

B5 = B4 + ones(7,1)*A(1,:);

%write_xyz(aa,[A;B5],'AB5.xyz');

% now rotate by 
