from subprocess import call

barcode = {
    'BAC01.BanII.cleaned.paired.1' : 'TGCGATAT',
    'BAC01.BanII.cleaned.paired.2' : 'CCTAT',
    'BAC02.BanII.cleaned.paired.1' : 'TTAGAGTGTG',
    'BAC02.BanII.cleaned.paired.2' : 'CCTAT',
    'BAC03.BanII.cleaned.paired.1' : 'ACGAGTG',
    'BAC03.BanII.cleaned.paired.2' : 'AGAGAT',
    'BAC04.BanII.cleaned.paired.1' : 'GGCTGTGTAT',
    'BAC04.BanII.cleaned.paired.2' : 'AGAGAT',
    'BAC05.BanII.cleaned.paired.1' : 'TATCTCGT',
    'BAC05.BanII.cleaned.paired.2' : 'TTGCGAT',
    'BAC06.BanII.cleaned.paired.1' : 'ACATACTC',
    'BAC06.BanII.cleaned.paired.2' : 'TTGCGAT',
    'BAC07.BanII.cleaned.paired.1' : 'CTGTGTATAC',
    'BAC07.BanII.cleaned.paired.2' : 'AACGTCTC',
    'BAC08.BanII.cleaned.paired.1' : 'GGTCTCAGT',
    'BAC08.BanII.cleaned.paired.2' : 'AACGTCTC',
    'BAC09.BanII.cleaned.paired.1' : 'CTACTCG',
    'BAC09.BanII.cleaned.paired.2' : 'GTACACCTA',
    'BAC10.BanII.cleaned.paired.1' : 'ACGATGCAT',
    'BAC10.BanII.cleaned.paired.2' : 'GTACACCTA',
    'BAC11.BanII.cleaned.paired.1' : 'CCTAGAC',
    'BAC11.BanII.cleaned.paired.2' : 'TACTCGACTG',
    'BAC12.BanII.cleaned.paired.1' : 'TGATGTGA',
    'BAC12.BanII.cleaned.paired.2' : 'TACTCGACTG',
    'BAC01.NlaIII.cleaned.paired.1' : 'GGATACTGT',
    'BAC01.NlaIII.cleaned.paired.2' : 'CCTAT',
    'BAC02.NlaIII.cleaned.paired.1' : 'CTCTACATGT',
    'BAC02.NlaIII.cleaned.paired.2' : 'CCTAT',
    'BAC03.NlaIII.cleaned.paired.1' : 'TTGTGTCGC',
    'BAC03.NlaIII.cleaned.paired.2' : 'AGAGAT',
    'BAC04.NlaIII.cleaned.paired.1' : 'TCTGCGT',
    'BAC04.NlaIII.cleaned.paired.2' : 'AGAGAT',
    'BAC05.NlaIII.cleaned.paired.1' : 'CATCATCAGC',
    'BAC05.NlaIII.cleaned.paired.2' : 'TTGCGAT',
    'BAC06.NlaIII.cleaned.paired.1' : 'CTACGAGT',
    'BAC06.NlaIII.cleaned.paired.2' : 'TTGCGAT',
    'BAC07.NlaIII.cleaned.paired.1' : 'TACTACGCT',
    'BAC07.NlaIII.cleaned.paired.2' : 'AACGTCTC',
    'BAC08.NlaIII.cleaned.paired.1' : 'ATCACGCAT',
    'BAC08.NlaIII.cleaned.paired.2' : 'AACGTCTC',
    'BAC09.NlaIII.cleaned.paired.1' : 'TCAGACTG',
    'BAC09.NlaIII.cleaned.paired.2' : 'GTACACCTA',
    'BAC10.NlaIII.cleaned.paired.1' : 'CGCACAT',
    'BAC10.NlaIII.cleaned.paired.2' : 'GTACACCTA',
    'BAC11.NlaIII.cleaned.paired.1' : 'TTGTGCAT',
    'BAC11.NlaIII.cleaned.paired.2' : 'TACTCGACTG',
    'BAC12.NlaIII.cleaned.paired.1' : 'GAGCATAGT',
    'BAC12.NlaIII.cleaned.paired.2' : 'TACTCGACTG',
    'BAC01.NspI.cleaned.paired.1' : 'TGCATGTCG',
    'BAC01.NspI.cleaned.paired.2' : 'CCTAT',
    'BAC02.NspI.cleaned.paired.1' : 'CGCATGTCAT',
    'BAC02.NspI.cleaned.paired.2' : 'CCTAT',
    'BAC03.NspI.cleaned.paired.1' : 'ATATGCTAGC',
    'BAC03.NspI.cleaned.paired.2' : 'AGAGAT',
    'BAC04.NspI.cleaned.paired.1' : 'TCTGTGAGC',
    'BAC04.NspI.cleaned.paired.2' : 'AGAGAT',
    'BAC05.NspI.cleaned.paired.1' : 'ATGTGCGATG',
    'BAC05.NspI.cleaned.paired.2' : 'TTGCGAT',
    'BAC06.NspI.cleaned.paired.1' : 'TCTGTCATGT',
    'BAC06.NspI.cleaned.paired.2' : 'TTGCGAT',
    'BAC07.NspI.cleaned.paired.1' : 'TCGTAGA',
    'BAC07.NspI.cleaned.paired.2' : 'AACGTCTC',
    'BAC08.NspI.cleaned.paired.1' : 'TCTAGAGA',
    'BAC08.NspI.cleaned.paired.2' : 'AACGTCTC',
    'BAC09.NspI.cleaned.paired.1' : 'TCACGCTGTA',
    'BAC09.NspI.cleaned.paired.2' : 'GTACACCTA',
    'BAC10.NspI.cleaned.paired.1' : 'TTAGTGCTC',
    'BAC10.NspI.cleaned.paired.2' : 'GTACACCTA',
    'BAC11.NspI.cleaned.paired.1' : 'CACGTCGTAT',
    'BAC11.NspI.cleaned.paired.2' : 'TACTCGACTG',
    'BAC12.NspI.cleaned.paired.1' : 'ACAGTCGT',
    'BAC12.NspI.cleaned.paired.2' : 'TACTCGACTG',
    'BAC01.Bsp1286Iset2.cleaned.paired.1' : 'CTAGTACAG',
    'BAC01.Bsp1286Iset2.cleaned.paired.2' : 'TTGCGAT',
    'BAC02.Bsp1286Iset2.cleaned.paired.1' : 'TCATAGACGC',
    'BAC02.Bsp1286Iset2.cleaned.paired.2' : 'TTGCGAT',
    'BAC03.Bsp1286Iset2.cleaned.paired.1' : 'TGTGATG',
    'BAC03.Bsp1286Iset2.cleaned.paired.2' : 'TTGCGAT',
    'BAC04.Bsp1286Iset2.cleaned.paired.1' : 'CGCTACTAGT',
    'BAC04.Bsp1286Iset2.cleaned.paired.2' : 'TTGCGAT',
    'BAC05.Bsp1286Iset2.cleaned.paired.1' : 'TGTGCATCAC',
    'BAC05.Bsp1286Iset2.cleaned.paired.2' : 'TTGCGAT',
    'BAC06.Bsp1286Iset2.cleaned.paired.1' : 'TCACTATC',
    'BAC06.Bsp1286Iset2.cleaned.paired.2' : 'TTGCGAT',
    'BAC07.Bsp1286Iset2.cleaned.paired.1' : 'CTAGTACAG',
    'BAC07.Bsp1286Iset2.cleaned.paired.2' : 'AACGTCTC',
    'BAC08.Bsp1286Iset2.cleaned.paired.2' : 'AACGTCTC',
    'BAC08.Bsp1286Iset2.cleaned.paired.1' : 'TCATAGACGC',
    'BAC09.Bsp1286Iset2.cleaned.paired.2' : 'AACGTCTC',
    'BAC09.Bsp1286Iset2.cleaned.paired.1' : 'TGTGATG',
    'BAC10.Bsp1286Iset2.cleaned.paired.2' : 'AACGTCTC',
    'BAC10.Bsp1286Iset2.cleaned.paired.1' : 'CGCTACTAGT',
    'BAC11.Bsp1286Iset2.cleaned.paired.1' : 'TGTGCATCAC',
    'BAC11.Bsp1286Iset2.cleaned.paired.2' : 'AACGTCTC',
    'BAC12.Bsp1286Iset2.cleaned.paired.1' : 'TCACTATC',
    'BAC12.Bsp1286Iset2.cleaned.paired.2' : 'AACGTCTC',
}

#RE = ["BanII"]
RE = ["BanII", "Bsp1286Iset2", "NspI", "NlaIII"]
#BAC = ["BAC01"]
BAC = ["BAC01","BAC02","BAC03","BAC04","BAC05","BAC06","BAC07","BAC08","BAC09","BAC10","BAC11","BAC12"]

for re in RE:
    for bac in BAC:
        r_1 = barcode[bac+"."+re+".cleaned.paired.1"][::-1]
        r_2 = barcode[bac+"."+re+".cleaned.paired.2"][::-1]
        r_1 = r_1.translate(str.maketrans('ATCGatcg','TAGCtagc'))
        r_2 = r_2.translate(str.maketrans('ATCGatcg','TAGCtagc'))

        call("cutadapt --info-file error1 -f fastq -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -g ^"+barcode[bac+"."+re+".cleaned.paired.1"]+" -g ^"+barcode[bac+"."+re+".cleaned.paired.2"]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+bac+"."+re+".tmp.1.fastq -o "+bac+"."+re+".tmp.2.fastq ../../input/re/"+re+"_B73_"+bac+"_"+barcode[bac+"."+re+".cleaned.paired.1"]+"-"+barcode[bac+"."+re+".cleaned.paired.2"]+".trimmed-paired-1.fq ../../input/re/"+re+"_B73_"+bac+"_"+barcode[bac+"."+re+".cleaned.paired.1"]+"-"+barcode[bac+"."+re+".cleaned.paired.2"]+".trimmed-paired-2.fq",shell=True)
        call("cutadapt --info-file error2 -f fastq -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^"+barcode[bac+"."+re+".cleaned.paired.1"]+" -g ^"+barcode[bac+"."+re+".cleaned.paired.2"]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+bac+"."+re+".cleaned.paired.1.cutadapter.fastq -o "+bac+"."+re+".cleaned.paired.2.cutadapter.fastq "+bac+"."+re+".tmp.1.fastq "+bac+"."+re+".tmp.2.fastq",shell=True)
        call("python3 head_trans_tail.py "+bac+"."+re+".cleaned.paired.1.cutadapter.fastq > "+bac+"."+re+".cleaned.paired.1.cutadapter.R.fastq",shell=True)
        call("python3 head_trans_tail.py "+bac+"."+re+".cleaned.paired.2.cutadapter.fastq > "+bac+"."+re+".cleaned.paired.2.cutadapter.R.fastq",shell=True)
        call("cutadapt --info-file error3 -f fastq -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -g ^"+barcode[bac+"."+re+".cleaned.paired.1"]+" -g ^"+barcode[bac+"."+re+".cleaned.paired.2"]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+bac+"."+re+".tmpagain.1.fastq -o "+bac+"."+re+".tmpagain.2.fastq "+bac+"."+re+".cleaned.paired.1.cutadapter.R.fastq "+bac+"."+re+".cleaned.paired.2.cutadapter.R.fastq",shell=True)
        call("cutadapt --info-file error4 -f fastq -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^"+barcode[bac+"."+re+".cleaned.paired.1"]+" -g ^"+barcode[bac+"."+re+".cleaned.paired.2"]+" -g ^"+r_1+" -g ^"+r_2+" -n 2 -e 0.15 --paired-output "+bac+"."+re+".cleaned.paired.1.cutadapter.2ends.fastq -o "+bac+"."+re+".cleaned.paired.2.cutadapter.2ends.fastq "+bac+"."+re+".tmpagain.1.fastq "+bac+"."+re+".tmpagain.2.fastq",shell=True)
        call("python3 head_trans_tail.py "+bac+"."+re+".cleaned.paired.1.cutadapter.2ends.fastq > "+bac+"."+re+".1.cutadapter.2ends.fastq",shell=True)
        call("python3 head_trans_tail.py "+bac+"."+re+".cleaned.paired.2.cutadapter.2ends.fastq > "+bac+"."+re+".2.cutadapter.2ends.fastq",shell=True)
        #call("mv $bac.$re.cleaned.paired.cutadapter.2ends.fastq $bac.$re.cutadapter.2ends.fastq",shell=True)
        #call("rm "+bac+"."+re+".cleaned.paired.*.cutadapter.* "+bac+"."+re+".tmp*.?.fastq",shell=True)
        call("mv "+bac+"."+re+".1.cutadapter.2ends.fastq output",shell=True)
        call("mv "+bac+"."+re+".2.cutadapter.2ends.fastq output",shell=True)
	#error correction#
	#call(["./musket -p 24 $bac.$re.cutadapter.2ends.fastq -omulti $bac.$re.cutadapter.2ends.ec.fastq -inorder"],shell=True)
	#call(["mv "+bac+"."+re+".cutadapter.2ends.ec.fastq.0 $bac.$re.cutadapter.2ends.ec.fastq")
        #call(["cutadapt -f fastq -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC"],shell=True)
