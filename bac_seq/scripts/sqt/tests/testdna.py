from sqt.dna import reverse_complement

def test_complement_string():
	rc = reverse_complement
	assert rc('') == ''
	assert rc('A') == 'T'
	assert rc('C') == 'G'
	assert rc('TG') == 'CA'
	assert rc('N') == 'N'
	assert rc('a') == 't'

	assert rc('ACGTUMRWSYKVHDBN') == 'NVHDBMRSWYKAACGT'
	assert rc('acgtumrwsykvhdbn') == 'nvhdbmrswykaacgt'
	assert rc('ACGTUMRWSYKVHDBNacgtumrwsykvhdbn') == 'nvhdbmrswykaacgtNVHDBMRSWYKAACGT'

#ACGTUMRWSYKVHDBN
#TGCAAKYWSRMBDHVN


def test_complement_bytes():
	rc = reverse_complement
	assert rc(b'') == b''
	assert rc(b'A') == b'T'
	assert rc(b'C') == b'G'
	assert rc(b'TG') == b'CA'
	assert rc(b'N') == b'N'
	assert rc(b'a') == b't'

	assert rc(b'ACGTUMRWSYKVHDBN') == b'NVHDBMRSWYKAACGT'
	assert rc(b'acgtumrwsykvhdbn') == b'nvhdbmrswykaacgt'
	assert rc(b'ACGTUMRWSYKVHDBNacgtumrwsykvhdbn') == b'nvhdbmrswykaacgtNVHDBMRSWYKAACGT'
