from sqt.align import edit_distance as ed, GlobalAlignment as GA

STRING_PAIRS = [
	(b'', b''),
	(b'A', b'A'),
	(b'A', b''),
	(b'AB', b'ABC'),
	(b'ANANAS', b'BANANA'),
	(b'SISSI', b'MISSISSIPPI'),
	(b'GGAATCCC', b'TGAGGGATAAATATTTAGAATTTAGTAGTAGTGTT'),
	(b'TCTGTTCCCTCCCTGTCTCA', b'TTTTAGGAAATACGCC'),
	(b'TGAGACACGCAACATGGGAAAGGCAAGGCACACAGGGGATAGG', b'AATTTATTTTATTGTGATTTTTTGGAGGTTTGGAAGCCACTAAGCTATACTGAGACACGCAACAGGGGAAAGGCAAGGCACA'),
	(b'TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA', b'TTTTAGGAAATACGCCTGGTGGGGTTTGGAGTATAGTGAAAGATAGGTGAGTTGGTCGGGTG'),
	(b'A', b'TCTGCTCCTGGCCCATGATCGTATAACTTTCAAATTT'),
	]


def test_edit_distance():
	assert ed(b'', b'') == 0
	assert ed(b'', b'A') == 1
	assert ed(b'A', b'B') == 1
	assert ed(b'A', b'A') == 0
	assert ed(b'A', b'AB') == 1
	assert ed(b'BA', b'AB') == 2
	for s, t in STRING_PAIRS:
		assert ed(s, b'') == len(s)
		assert ed(b'', s) == len(s)
		assert ed(s, t) == ed(t, s)


def nongap_characters(row):
	"""
	Return the non-gap characters (not '\0') of an alignment row.
	"""
	try:
		return row.replace(b'\0', b'')
	except TypeError:
		return row.replace('\0', '')


def count_gaps(row):
	try:
		return row.count(b'\0')
	except TypeError:
		return row.count('\0')


def count_mismatches(row1, row2):
	if type(row1) is str:
		gap = '\0'
	else:
		gap = 0
	return sum(1 for (c1, c2) in zip(row1, row2) if c1 != c2 and c1 != gap and c2 != gap)


def test_global_alignment():
	for s, t in STRING_PAIRS:
		distance = ed(s, t)
		ga = GA(s, t)
		assert len(ga.row1) == len(ga.row2)
		assert ga.errors == distance
		assert nongap_characters(ga.row1) == s
		assert nongap_characters(ga.row2) == t
		assert ga.errors == count_gaps(ga.row1) + count_gaps(ga.row2) + count_mismatches(ga.row1, ga.row2)

