#!/bin/bash

for B in {BAC01,BAC02,BAC03,BAC04,BAC05,BAC06,BAC07,BAC08,BAC09,BAC10,BAC11,BAC12}
do
    python readcut_NlaIII.py ../output/$B.BanII.1.cutadapter.2ends.fastq ../output/$B.BanII.2.cutadapter.2ends.fastq > check2
    echo $(wc -l check2 | awk '{print $1}' | head -1)/9 | bc
done
