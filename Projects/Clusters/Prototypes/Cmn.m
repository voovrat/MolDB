function cmn=Cmn(m,n)

if m==0 || m==n
   cmn = 1;
   return
end

cmn=Cmn(m,n-1) + Cmn(m-1,n-1);
