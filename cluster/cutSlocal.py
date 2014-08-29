#run as <python ../../bowtieindex/sam_file>

import sys
import re
i=0
f1 = open(sys.argv[1],'r')
lines = f1.readlines()
f1.close()

for line in lines:
    frontS=0
    endS=0
    if i<4:
        i=i+1
    else:
        parts = line.split()
        valid = parts[5]
        if valid!='*':
            print "@"+parts[0]
            a=re.split("S*\d+[DIM]",valid)
            if a[0]!='':
                frontS = int(a[0])
            if a[-1]!='':
                endS = int(a[-1][0:-1])

            if endS!=0:
                print parts[9][frontS:-endS]+"\n+\n"+parts[10][frontS:-endS]
            else:
                print parts[9][frontS:]+"\n+\n"+parts[10][frontS:]
