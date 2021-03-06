[A,at] = read_xyz('SPCE.xyz');

D = mydist_fast(A);

oh = ( D(1,2) + D(1,3) )/2;
hh = D(2,3);

%  hh^2 = 2oh^2 -2oh^2 cos th ==>
% 1 - cos th = 0.5 *(hh/oh)^2
% cos th = 1 - 0.5 * (hh/oh)^2 

cos_th = 1 - 0.5 * (hh/oh)^2
th = acos(cos_th)



th*180/pi % 109


B = [ 0    0        0;
      0  sin(th/2) cos(th/2);
      0 -sin(th/2) cos(th/2) ]

write_xyz(at,B,'SPCE_sym.xyz')

D1 = mydist_fast(B)

D - D1;

   % ------------------------------ %



a = 5;   % 0.8 0.9 1 1.1 1.3 1.5 1.7 2 2.5 3 --> 10
b = 2.2;   % 0:8 0.9 1 1.1 1.3 1.5 1.7 2 2.5 3 --> 10
c = 0.8;   % 0 0.1 0.2 0.4 0.6 0.8 1 1.5 2 2.5 --> 10

th1 = 0;  % 0 pi/4 pi/2 3pi/4 pi         -->5 
phi1 = pi/2; % 0 pi/4 pi/2 3pi/4 pi 5pi/4 3pi/2 7pi/4 --> 8
psi1 = 0;  % 0 pi/4 pi/2 3pi/4   --> 4 
% except for th=0,pi one have always phi=0

%% 5 * 8 * 4 = 160
%% in clever case (account for specials for th=0,pi:
%  3*32 + 2*4 --> 104


th2 = pi;
phi2= 0;
psi2 = 0;


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


[W0,at] = read_xyz('SPCE_sym.xyz');

W1 = (Rot1 * W0')' + ones(3,1) * [0 0 a];
W2 = (Rot2 * W0')' - ones(3,1) * [0 0 b];

p = [ c 0 0];

at2 = [at;at; {'H'} ];

write_xyz(at2,[W1;W2;p],'sys.xyz');










