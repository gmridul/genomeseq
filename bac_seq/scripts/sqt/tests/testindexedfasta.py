from nose.tools import raises
from sqt.io.fasta import IndexedFasta, FastaReader
import os.path

def dpath(path):
	return os.path.join(os.path.dirname(__file__), path)


@raises(ValueError)
def test_indexedfasta_contextmanager():
	indfasta = IndexedFasta(dpath("seq.fa"))
	with indfasta as ifw:
		pass
	with indfasta as ifw:
		pass


def test_indexedfasta():
	with IndexedFasta(dpath("seq.fa")) as ifa:
		assert len(ifa) == 3
		chr1 = ifa.get("Chr1")
		chr2 = ifa.get("Chr2")
		assert len(chr1) == 1235
		assert chr1[0:300].startswith(b'CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATC')
		assert chr1[:300].startswith(b'CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATC')
		assert chr1[:].startswith(b'CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATC')
		assert chr2[227:320] == b'gttggaatcgTTCCGAGTTTTCTCAGCAGTTCTCGGACAAAAACTGATGAATCGTCGAGGAGAATGAGCTTGCCTTGCGTGGGCTGCCATTAG'
		assert chr1[:300].startswith(b'CCCTAAACCCTA')
		assert chr2[:].endswith(b'TATCCGAGGGATGGTATCGG')


def test_all_regions():
	# read the file via a FastaReader, then check that all substrings are equal
	path = dpath("indexed.fasta")
	sequences = dict()
	with FastaReader(path, binary=True) as fr:
		for record in fr:
			sequences[record.name] = record.sequence
	with IndexedFasta(path):
		indexed = IndexedFasta(path)
	
	regions = []
	for name in sorted(sequences):
		for i in range(len(sequences[name])):
			for j in range(i, len(sequences[name])):
				regions.append( (name, i, j) )

	for region in regions:
		expected = sequences[region[0]][region[1]:region[2]]
		i = indexed.get(region[0])
		print("i:", i)
		assert indexed.get(region[0])[region[1]:region[2]] == expected
