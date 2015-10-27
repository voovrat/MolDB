function build_water(l1,l2,alpha,output)



XYZ = [ 0 0 0; ...
        l1*cos(alpha/2)  l1*sin(alpha/2) 0; ...
       -l2*cos(alpha/2)  l2*sin(alpha/2) 0 ];

at = {'O';'H';'H'};

write_xyz(at,XYZ,output);


