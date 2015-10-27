function build_water_pair(d,outfile)
%  s = struct( a,b,c, th1, phi1, psi1, th2, phi2,psi2 )
%



a = d.a;   % 0.8 0.9 1 1.1 1.3 1.5 1.7 2 2.5 3 --> 10
b = d.b;   % 0:8 0.9 1 1.1 1.3 1.5 1.7 2 2.5 3 --> 10
c = d.c;   % 0 0.1 0.2 0.4 0.6 0.8 1 1.5 2 2.5 --> 10

th1 = d.th1;  % 0 pi/4 pi/2 3pi/4 pi         -->5 
phi1 = d.phi1; % 0 pi/4 pi/2 3pi/4 pi 5pi/4 3pi/2 7pi/4 --> 8
psi1 = d.psi1;  % 0 pi/4 pi/2 3pi/4   --> 4 
% except for th=0,pi one have always phi=0

%% 5 * 8 * 4 = 160
%% in clever case (account for specials for th=0,pi:
%  3*32 + 2*4 --> 104


th2 = d.th2;
phi2= d.phi2;
psi2 = d.psi2;


%% --> 10*10*10*104^2 = 10.8 M

% rotate order: psi Oz - theta Oy - phi Oz

RotPsi1 = [ cos(psi1) -sin(psi1) 0;
            sin(psi1)  cos(psi1) 0;
               0        0      1 ];

RotTh1 = [ cos(th1)    0   -sin(th1) ;
              0        1      0  ;
           sin(th1)    0   cos(th1) ];

RotPhi1 =  [ cos(phi1) -sin(phi1) 0;
             sin(phi1)  cos(phi1) 0;
                 0          0     1];

Rot1= RotPhi1 * RotTh1 * RotPsi1;


%%
RotPsi2 = [ cos(psi2) -sin(psi2) 0;
            sin(psi2)  cos(psi2) 0;
               0        0      1 ];

RotTh2 = [ cos(th2)    0   -sin(th2) ;
              0        1      0  ;
           sin(th2)    0   cos(th2) ];

RotPhi2 =  [ cos(phi2) -sin(phi2) 0;
             sin(phi2)  cos(phi2) 0;
                 0          0     1];

Rot2= RotPhi2 * RotTh2 * RotPsi2;


W0 = [0.0000000000 0.0000000000 0.0000000000
 0.0000000000 0.8141155184 0.5807029557
 0.0000000000 -0.8141155184 0.5807029557 ];

at = {'O';'H';'H'};


%[W0,at] = read_xyz('SPCE_sym.xyz');

W1 = (Rot1 * W0')' + ones(3,1) * [0 0 a];
W2 = (Rot2 * W0')' - ones(3,1) * [0 0 b];

p = [ c 0 0];

at2 = [at;at; {'H'} ];

write_xyz(at2,[W1;W2;p],outfile);



