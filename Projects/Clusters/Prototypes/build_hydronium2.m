function build_hydronium2(l1,l2,l3,alpha,beta,gamma,output);
%        x 
%       H1                         /H3
%       |                 _____l3/_|___x
% alpha | l1             /   H1-O gamma
% y-    O               /     / \ | /
%   l2/  \ l3cos(gamma)/____H2__\|_/ y
%    /beta\                 H3proj
%   H2     H3proj
%
% alpha = 0..pi (around 120 deg)
% beta = 0..pi  (around 120 deg )
% gamma = -pi/2..pi/2 (around 0 )
%
% l1,l2,l3 = 0..inf ( around 1.1 )
%

H1 = [  l1 0  0]
H2 = [ cos(alpha)*l2  sin(alpha)*l2 0]
H3proj = [ cos(alpha+beta) * l3*cos(gamma) sin(alpha+beta) * l3*cos(gamma) 0]

H3 = [ cos(alpha+beta) * l3*cos(gamma) ...
       sin(alpha+beta)*l3*cos(gamma) ...
       sin(gamma)*l3 ]

at = {'O','H','H','H'}';
XYZ = [ 0 0 0; H1;H2;H3];

write_xyz(at,XYZ,output);


