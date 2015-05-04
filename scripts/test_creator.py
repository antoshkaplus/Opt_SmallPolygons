

import sys
import os

N = int(sys.argv[1])
for i in range(1, N+1):
    os.system("java -jar tester.jar -exec 'python test_writer.py ../data/output_%d.txt' -novis -seed %d" % (i, i))
    
    
    