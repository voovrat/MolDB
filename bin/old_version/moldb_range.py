import sys,os


if __name__ == '__main__':

    if len(sys.argv)<3:
        print "Useage:  moldb_range   from to step"
        quit();

    From = float(sys.argv[1])
    To = float(sys.argv[2]);
    Step = float(sys.argv[3]);

    eps = 1e-8

    x = From;
    while x<To+eps:
        print  "%g"%x
        x+=Step
