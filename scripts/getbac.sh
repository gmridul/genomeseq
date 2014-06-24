for i in AC217266 #AC217290 #AC217267 AC217268 AC217269 AC217270 AC217271 AC217272 AC217273 AC217274 AC217275 AC217276 AC217277 AC217278 AC217279 AC217280 AC217281 AC217282 AC217283 AC217284 AC217285 AC217286 AC217287 AC217288 AC217289 
do
    wget http://www.ncbi.nlm.nih.gov/nuccore/$i
    grep ReportShortCut15 $i > line$i 
    rm $i
    echo http://www.ncbi.nlm.nih.gov$(awk '{print $3}' line$i | awk 'BEGIN { FS="\"" } {print $2}' | head -1)\&format=text
    casperjs crawl.js http://www.ncbi.nlm.nih.gov$(awk '{print $3}' line$i | awk 'BEGIN { FS="\"" } {print $2}' | head -1)\&format=text > $i
    #rm line$i
    sed -i -e 1,1d $i
done
