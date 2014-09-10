#copy the script on command line and replace 68443 with the required number ( it is the position at which the reads starts aligning to the reference genome"
awk  '{
if ($4 == 36844 && and($2,4)==0) {
    printf "@"; print $1; 
    if(and($2,16)==16) {
        len = length($10)
        for(i=len;i!=0;i--) {
            if(substr($10,i,1)=="'"A"'") printf "'"T"'"
            else if(substr($10,i,1)=="'"T"'") printf "'"A"'"
            else if(substr($10,i,1)=="'"C"'") printf "'"G"'"
            else if(substr($10,i,1)=="'"G"'") printf "'"C"'"
        }
        print "\n+"; 
        x=""
        for(i=len;i!=0;i--) {
            printf substr($11,i,1)
        }
        print ""
    }
    else {
        print $10; print "+"; print $11
    }
}}' < $1
