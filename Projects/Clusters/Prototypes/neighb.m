
function X = neighb(X0,R)

   
   D = sqrt( sum(X0(1:3:end,:).^2,2));
   I0 = find(D<R);
   [DD,ID] = sort(D(I0));
   I = I0(ID);
%   II = [ 3*I - 2; 3*I-1; 3*I];
  Nox = size(I,1);
  II=zeros(3*Nox,1);
  II(1:3:end) = 3*I-2;
  II(2:3:end) = 3*I-1;
  II(3:3:end) = 3*I;
   X = X0(II,:);

