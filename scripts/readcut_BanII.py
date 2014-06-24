import sys
import re
from itertools import izip

a=re.compile("G[AG]GC[CT]C")
with open(sys.argv[1],'r') as f1:
    with open(sys.argv[2],'r') as f2:
        for l in izip(f1,f1,f1,f1,f2,f2,f2,f2):
            if a.search(l[1]) or a.search(l[5]):
                print ''.join(l)
