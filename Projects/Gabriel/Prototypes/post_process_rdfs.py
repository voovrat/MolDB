import sys

def transp(matrix):
  n = len(matrix[0])
  m = len(matrix)
  tlines = [ [] for j in range(n) ]
  for i in range(m):
    for j in range(n):
      tlines[j].append(matrix[i][j])
  return tlines
  
def getMatrix(lines):
  m=len(lines)
  matrix=[]
  for l in lines:
    matrix.append( [ float(w) for w in l.split() ])
  return matrix


if __name__ == '__main__':

  if len(sys.argv)<=3:
    print "Usage: python post_process_rdfs.py  file.gro0  file.rdfs output_prefix"

  gro0_file = sys.argv[1]
  rdfs_file = sys.argv[2]
  prefix = sys.argv[3]

  f = open(gro0_file)
  gro0_lines = f.readlines()
  f.close()

  f = open(rdfs_file)
  rdfs_lines = f.readlines()
  f.close()


  rdfs = transp(getMatrix(rdfs_lines))

#  print rdfs

  atom_index={}

  for i in range(len(gro0_lines)):
    l=gro0_lines[i]
    w = l.split()
    at=w[1]
    if not at in atom_index: 
      atom_index[at] = [ i+1 ]  # i+1 because 0th column is r
    else:
      atom_index[at].append(i+1)

  atom_rdfs = {}

  for at in atom_index:
    idx = atom_index[at]
    NR = len(idx)
    N = len(rdfs[0])
    SumG = [ 0 for i in range(N) ]
    for i in idx:
      for j in range(N):
        SumG[j] += rdfs[i][j]
    atom_rdfs[at] = [ gi/NR for gi in SumG ]

  for at in atom_rdfs:
    f=open(prefix + at+'.rdf','w');
    for i in range(N):
      f.write(str(rdfs[0][i]) + ' '+str(atom_rdfs[at][i])+'\n')
    f.close()

  
