from sqt.cigar import parse, Cigar

def test_parse():
	assert parse("4S17M8D4M9I3H") == [(4, 4), (0, 17), (2, 8), (0, 4), (1, 9), (5, 3)]


def test_cigar_class():
	assert Cigar('4M') == Cigar([(0, 4)])
	c = Cigar('4S 17M 8D 4M 9I 3H')
	assert str(c) == '4S17M8D4M9I3H'
	assert '{}'.format(c) == str(c)
	assert '{: }'.format(c) == '4S 17M 8D 4M 9I 3H'
	assert Cigar('4M') + Cigar('1D') == Cigar('4M 1D')
	assert Cigar('2S 4M') + Cigar('3M') == Cigar('2S 7M')

