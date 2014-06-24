from sqt.scripts.fastaextract import parse_region

def test_parse_region():
	assert parse_region("chr7:5-7") == ("chr7", 4, 7, False)
	assert parse_region("chr7") == ("chr7", 0, None, False)
	assert parse_region("rc:chr7") == ("chr7", 0, None, True)
	assert parse_region("chr7:1-100") == ("chr7", 0, 100, False)
	assert parse_region("rc:chr7:1-100") == ("chr7", 0, 100, True)
	assert parse_region("chr7:1-") == ("chr7", 0, None, False)
	assert parse_region("rc:chr7:1-") == ("chr7", 0, None, True)
