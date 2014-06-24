import sys
from itertools import izip

with open (sys.argv[1],'r') as f1:
    with open(sys.argv[2],'r') as f2:
        for l in izip(f1,f1,f2,f2,f2,f2):
            if not l[1].startswith(l[3]):
                print l[1],l[3]
