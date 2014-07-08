import sys
from itertools import izip
count=0
with open(sys.argv[1],'r') as f1:
        for l in izip(f1,f1,f1,f1):
            count=count+1
            if len(l[1]) != len(l[3]):
                print(count)
