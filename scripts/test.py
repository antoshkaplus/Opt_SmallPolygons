

import sys
import os

N = int(sys.argv[1])
for i in range(1, N+1):
    print "seed", i
    os.system("java -jar tester.jar -exec 'python test_runner.py ../data/in_test_%d.txt ../data/out_test_%d.txt' -novis -seed %d" % (i, i, i))
    print
    print

