"""
Tests for the sqt.io.fasta module
"""
from nose.tools import raises

from sqt.io.fasta import FastaReader, FastaWriter, Sequence, FastqWriter
import os.path


def dpath(path):
	return os.path.join(os.path.dirname(__file__), path)


def test_fastqwriter():
	tmp = dpath("tmp.fastq")
	with FastqWriter(tmp) as fq:
		fq.write("name", "CCATA", "!#!#!")
		fq.write("name2", "HELLO", "&&&!&&")
	assert fq._file.closed
	with open(tmp) as t:
		assert t.read() == '@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n'
	os.remove(tmp)

def test_fastqwriter_twoheaders():
	tmp = dpath("tmp.fastq")
	with FastqWriter(tmp, twoheaders=True) as fq:
		fq.write("name", "CCATA", "!#!#!")
		fq.write("name2", "HELLO", "&&&!&&")
	assert fq._file.closed
	with open(tmp) as t:
		assert t.read() == '@name\nCCATA\n+name\n!#!#!\n@name2\nHELLO\n+name2\n&&&!&&\n'
	os.remove(tmp)


def test_fastawriter():
	tmp = dpath("tmp.fasta")
	with FastaWriter(tmp) as fw:
		fw.write("name", "CCATA")
		fw.write("name2", "HELLO")
	assert fw._file.closed
	with open(tmp) as t:
		assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'
	os.remove(tmp)


def test_fastawriter_linelength():
	tmp = dpath("tmp.fasta")
	with FastaWriter(tmp, line_length=3) as fw:
		fw.write("name", "CCAT")
		fw.write("name2", "TACCAG")
	assert fw._file.closed
	with open(tmp) as t:
		d = t.read()
		assert d == '>name\nCCA\nT\n>name2\nTAC\nCAG\n'
	os.remove(tmp)


def test_fastawriter_sequence():
	tmp = dpath("tmp.fasta")
	with FastaWriter(tmp) as fw:
		fw.write(Sequence("name", "CCATA"))
		fw.write(Sequence("name2", "HELLO"))
	assert fw._file.closed
	with open(tmp) as t:
		assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'
	os.remove(tmp)


@raises(ValueError)
def test_fastawriter_contextmanager():
	tmp = dpath("tmp.fasta")
	fr = FastaWriter(tmp)
	os.remove(tmp)
	with fr as frw:
		pass
	with fr as frw:
		pass


def test_fastareader():
	with FastaReader(dpath("seq.fa")) as fr:
		seqs = list(fr)
	assert fr.fp.closed
	assert len(seqs) == 3
	assert seqs[0].qualities is None
	assert seqs[0].name == 'Chr1'
	assert seqs[1].name == 'Chr2 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02'
	assert len(seqs[0].sequence) == 1235
	assert seqs[0].sequence.startswith('CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATC')
	assert seqs[1].sequence.startswith('ctcgaccaggacgatgaatgggc')
	assert seqs[2].sequence.endswith('AATCTTGCAAGTTCCAACTAATT')

def test_fastareader_binary():
	with FastaReader(dpath("seq.fa"), binary=True) as fr:
		seqs = list(fr)
	assert fr.fp.closed
	assert len(seqs) == 3
	assert seqs[0].qualities is None
	assert seqs[0].name == 'Chr1'
	assert seqs[2].name == 'Chr3 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02'
	assert len(seqs[0].sequence) == 1235
	assert seqs[0].sequence.startswith(b'CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATC')
	assert seqs[1].sequence.startswith(b'ctcgaccaggacgatgaatgggc')
	assert seqs[2].sequence.endswith(b'AATCTTGCAAGTTCCAACTAATT')


@raises(ValueError)
def test_fastareader_contextmanager():
	fr = FastaReader(dpath("seq.fa"))
	with fr as frw:
		pass
	with fr as frw:
		pass
