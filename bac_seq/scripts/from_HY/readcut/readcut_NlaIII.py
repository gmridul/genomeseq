import sys
import re
from itertools import izip
count=0
a=re.compile("[NACTG]CATG[NACTG]")
with open(sys.argv[1],'r') as f1:
#    with open(sys.argv[2],'r') as f2:
        for l in izip(f1,f1,f1,f1):
            if (a.search(l[1])):
                count=count+1
 #               print ''.join(l)

print count
