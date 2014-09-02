#run as <python cutSlocal.py ../../bowtieindex/sam_file>

import sys
import re
i=0
f1 = open(sys.argv[1],'r') #pass the sam file as command line argument
lines = f1.readlines()
f1.close()

for line in lines:
    frontS=0    #stores the number of substitutions at the beginning of the read
    endS=0      #stores the number of substitutions at the end of the read
    if i<4:     #ignore the first 4 lines ( extra information in the input file )
        i=i+1
    else:
        parts = line.split()
        valid = parts[5]    #stores the information regarding mismatches (M,I,D,S)

#ignore if did match completely. * represents not matched
        if valid!='*':
            a=re.split("S*\d+[DIM]",valid)
            if a[0]!='':
                frontS = int(a[0])
            if a[-1]!='':
                endS = int(a[-1][0:-1])

            if frontS!=0:
                print "@"+parts[0]+":1"
                print parts[9][0:frontS]+"\n+\n"+parts[10][0:frontS]
            if endS!=0:
                print "@"+parts[0]+":2"
                print parts[9][-endS:]+"\n+\n"+parts[10][-endS:]
