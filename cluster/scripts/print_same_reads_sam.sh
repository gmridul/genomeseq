#copy the script on command line and replace 68443 with the required number ( it is the position at which the reads starts aligning to the reference genome"
awk  '{
if ($4 == 3854 && and($2,4)==0) {
    printf "@"; print $1; 
    #print "##########################################################"
    #print $2;
    #print $10;
    #print "##########################################################"
    if(and($2,16)==16) {
        len = length($10)
        print "**********************************************"
        for(i=len;i!=0;i--) {
            if(substr($10,i,1)=="'"A"'") printf "'"T"'"
            else if(substr($10,i,1)=="'"T"'") printf "'"A"'"
            else if(substr($10,i,1)=="'"C"'") printf "'"G"'"
            else if(substr($10,i,1)=="'"G"'") printf "'"C"'"
        }
        print "\n+"; 
        for(i=len;i!=0;i--) {
            printf substr($11,i,1)
        }
        print ""
    }
    else {
        print $10; print "+"; print $11
    }
}}' < $1
