from subprocess import call
import sys


barcode = {
 'AC185603.CATCATCAGC-CCTAT-1'    : 'CATCATCAGC',
 'AC186236.ACTGATC-CCTAT-1'       : 'ACTGATC',
 'AC189789.ATGTGCGATG-CCTAT-1'    : 'ATGTGCGATG',
 'AC194009.ACGCTAGT-CCTAT-1'      : 'ACGCTAGT',
 'AC205693.ACAGCATCTG-CCTAT-1'    : 'ACAGCATCTG',
 'AC205879.AGTCGAG-CCTAT-1'       : 'AGTCGAG',
 'AC206167.GCGTCAT-CCTAT-1'       : 'GCGTCAT',
 'AC208225.CAGATAGTGC-CCTAT-1'    : 'CAGATAGTGC',
 'AC209175.CACGTAC-CCTAT-1'       : 'CACGTAC',
 'AC213049.TTGTGCAT-CCTAT-1'      : 'TTGTGCAT',
 'AC214831.TACGTCGTG-CCTAT-1'     : 'TACGTCGTG',
 'AC217266.ACGAGTG-CCTAT-1'       : 'ACGAGTG',
}

barcode1 = {
 'AC185603.CATCATCAGC-CCTAT-1'    : 'CCTAT',
 'AC186236.ACTGATC-CCTAT-1'       : 'CCTAT',
 'AC189789.ATGTGCGATG-CCTAT-1'    : 'CCTAT',
 'AC194009.ACGCTAGT-CCTAT-1'      : 'CCTAT',
 'AC205693.ACAGCATCTG-CCTAT-1'    : 'CCTAT',
 'AC205879.AGTCGAG-CCTAT-1'       : 'CCTAT',
 'AC206167.GCGTCAT-CCTAT-1'       : 'CCTAT',
 'AC208225.CAGATAGTGC-CCTAT-1'    : 'CCTAT',
 'AC209175.CACGTAC-CCTAT-1'       : 'CCTAT',
 'AC213049.TTGTGCAT-CCTAT-1'      : 'CCTAT',
 'AC214831.TACGTCGTG-CCTAT-1'     : 'CCTAT',
 'AC217266.ACGAGTG-CCTAT-1'       : 'CCTAT',
}

for f in barcode:
    f2=f[:-1]+'2'
    r_1 = barcode[f][::-1]
    r_2 = barcode1[f][::-1]
    r_1 = r_1.translate(str.maketrans('ATCGatcg','TAGCtagc'))
    r_2 = r_2.translate(str.maketrans('ATCGatcg','TAGCtagc'))

    call("../../bac_seq/scripts/from_HY/cutadapt-1.4.2/bin/cutadapt --info-file error1 -f fastq -g ^AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -g ^"+barcode[f]+" -g ^"+barcode1[f]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+f+".tmp.fastq -o "+f2+".tmp.fastq ../"+f+".fq ../"+f2+".fq",shell=True)
    call("../../bac_seq/scripts/from_HY/cutadapt-1.4.2/bin/cutadapt --info-file error2 -f fastq -g ^AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^"+barcode[f]+" -g ^"+barcode1[f]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+f+".ncutadapter.fastq -o "+f2+".ncutadapter.fastq "+f+".tmp.fastq "+f2+".tmp.fastq",shell=True)
    call("python3 head_trans_tail.py "+f+".ncutadapter.fastq > "+f+".cutadapter.R.fastq",shell=True)
    call("python3 head_trans_tail.py "+f2+".ncutadapter.fastq > "+f2+".cutadapter.R.fastq",shell=True)
    call("../../bac_seq/scripts/from_HY/cutadapt-1.4.2/bin/cutadapt --info-file error3 -f fastq -g ^AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -g ^"+barcode[f]+" -g ^"+barcode1[f]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+f+".tmpagain.fastq -o "+f2+".tmpagain.fastq "+f+".cutadapter.R.fastq "+f2+".cutadapter.R.fastq",shell=True)
    call("../../bac_seq/scripts/from_HY/cutadapt-1.4.2/bin/cutadapt --info-file error4 -f fastq -g ^AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^"+barcode[f]+" -g ^"+barcode1[f]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+f+".cleaned.paired.cutadapter.2ends.fastq -o "+f2+".cleaned.paired.cutadapter.2ends.fastq "+f+".tmpagain.fastq "+f2+".tmpagain.fastq",shell=True)
    call("python3 head_trans_tail.py "+f+".cleaned.paired.cutadapter.2ends.fastq > "+f+".cutadapter.2ends.fastq",shell=True)
    call("python3 head_trans_tail.py "+f2+".cleaned.paired.cutadapter.2ends.fastq > "+f2+".cutadapter.2ends.fastq",shell=True)

    #call("mv $bac.$re.cleaned.paired.cutadapter.2ends.fastq $bac.$re.cutadapter.2ends.fastq",shell=True)
    call("rm *.cleaned.paired.* *.tmp*.* *.R.* *ncutadapter*",shell=True)
    #call("mv "+bac+"."+re+".1.cutadapter.2ends.fastq output",shell=True)
    #call("mv "+bac+"."+re+".2.cutadapter.2ends.fastq output",shell=True)
    #error correction#
    #call(["./musket -p 24 $bac.$re.cutadapter.2ends.fastq -omulti $bac.$re.cutadapter.2ends.ec.fastq -inorder"],shell=True)
    #call(["mv "+bac+"."+re+".cutadapter.2ends.ec.fastq.0 $bac.$re.cutadapter.2ends.ec.fastq")
    #call(["cutadapt -f fastq -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC"],shell=True)
