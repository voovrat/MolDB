function build_hydronium(z1,x2,z2,x3,y3,z3,output)
%
% build_hydronium(z1,x2,z2, x3,y3,z3,output) 
%
%           H1(0,0,z1)  
%           | 
%           |
%           O 
%          / \
%         /   H2(x2,0,z2)
%        /
%       H3(x3,y3,z3)
%


at = {'O';'H';'H';'H'};
XYZ = [ 0 0 0; 0 0 z1; x2 0 z2; x3 y3 z3];

write_xyz(at,XYZ,output);

