def combinations(m,n):
   if m==0:
      return [[]];
   if m==n: 
      return [range(n)];
   C = combinations(m,n-1)
   C1 = combinations(m-1,n-1)
   for c in C1:
      c.append(n-1)
   return C+C1

cc=[]

for i in range(1,8):
   c = combinations(i,7)
   cc += c

def cluster2str(c):
   s=''
   for i in range(len(c)):
      s+='c'+str(c[i])
   return s


for c in cc:
   print cluster2str(c)


