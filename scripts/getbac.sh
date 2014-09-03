cd ../input
for i in AC185603 AC186236 AC189789 AC194009 AC205693 AC205879 AC206167 AC208225 AC209175 AC213049 AC214831 AC217266
do
    wget http://www.ncbi.nlm.nih.gov/nuccore/$i 
    grep ReportShortCut15 $i > line$i 
    rm $i
    echo line$i
    echo http://www.ncbi.nlm.nih.gov$(cut -d"\"" -f4 line$i)\&format=text
    casperjs ../scripts/crawl.js http://www.ncbi.nlm.nih.gov$(cut -d"\"" -f4 line$i)\&format=text > $i
    rm line$i
    #sed -i -e 1,1d $i
done
