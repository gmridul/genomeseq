import sys
count=0
with open(sys.argv[1],'r') as fp:
    for line in fp:
        if count==1 or count==3:
           line=line[::-1]
           if count==1:
               line = line.translate(str.maketrans('ATCGatcg','TAGCtagc'))
        print(line.strip())
        count=(count+1)%4
