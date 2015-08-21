from getopt import gnu_getopt
import sys

if __name__=='__main__':
	optlist,args=gnu_getopt(sys.argv[1:],'abc')

	print optlist
	print args
