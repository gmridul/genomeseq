cd ../input
for i in AC217266
do
    wget http://www.ncbi.nlm.nih.gov/nuccore/$i 
    grep ReportShortCut15 $i > line$i 
    rm $i
    echo line$i
    echo http://www.ncbi.nlm.nih.gov$(awk '{print $3}' line$i | awk 'BEGIN { FS="\"" } {print $2}' | head -1)\&format=text
    ../scripts/casperjs crawl.js http://www.ncbi.nlm.nih.gov$(awk '{print $3}' line$i | awk 'BEGIN { FS="\"" } {print $2}' | head -1)\&format=text > $i
    rm line$i
    sed -i -e 1,1d $i
done
