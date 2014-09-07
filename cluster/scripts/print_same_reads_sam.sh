#copy the script on command line and replace 68443 with the required number ( it is the position at which the reads starts aligning to the reference genome"
awk  '{if ($4 == 68443) {printf "@"; print $1; print $10; print "+"; print $11}}' < AC217266.sam
