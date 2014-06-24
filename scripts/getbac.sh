for i in AC217265 AC217290 #AC217267 AC217268 AC217269 AC217270 AC217271 AC217272 AC217273 AC217274 AC217275 AC217276 AC217277 AC217278 AC217279 AC217280 AC217281 AC217282 AC217283 AC217284 AC217285 AC217286 AC217287 AC217288 AC217289 
do
    wget http://www.ncbi.nlm.nih.gov/nuccore/$i
    grep ReportShortCut15 $i > line$i 
    rm $i
    firefox http://www.ncbi.nlm.nih.gov$(awk '{print $3}' line$i | awk 'BEGIN { FS="\"" } {print $2}' | head -1)\&format=text &
    sleep 4
    xdotool key ctrl+a ctrl+c ctrl+w Return
    gedit $i &
    sleep 1
    xdotool key ctrl+v ctrl+s ctrl+w Return
    sleep 1
    rm line$i

done
