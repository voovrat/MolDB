from  math import sin,cos,pi
import sys


if __name__ == "__main__":
  
  if len(sys.argv)<3:
    print 'usage: python build_water l1 l2 alpha > output.xyz'
    print 'alpha in deg '
    quit()

  l1 = float( sys.argv[1] );
  l2 = float( sys.argv[2] );
  alpha = float( sys.argv[3] ) * pi / 180;

#  print 'l1='+str(l1) + ' l2='+str(l2) + ' alpha='+str(alpha)

  print '3'
  print 'water'
  print 'O 0.0 0.0 0.0'
  print 'H ' + str(  l1 * cos(pi/2 - alpha/2) ) + " " + str( l1 * sin(pi/2 - alpha/2) ) + " 0.0"
  print 'H ' + str( -l2 * cos(pi/2 - alpha/2) ) + " " + str( l2 * sin(pi/2 - alpha/2) ) + " 0.0"

# use 0.0, because Gaussian is too stupid to recognize that 0.0 and 0 is the same! 

