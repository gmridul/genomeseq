from sqt.math import frequency_median as median

def test_median():
	assert median( { 5: 2, 8: 4 } ) == 8
	assert median( { 5: 2, 8: 3 } ) == 8
	assert median( { 5: 1, 19: 2 } ) == 19
	assert median( { 27: 20, 5: 1, 19: 2 } ) == 27

	# one value
	assert median( { 5: 1 } ) == 5
	assert median( { 5: 1000 } ) == 5

	# five values
	assert median( { 5: 0, 8: 5 } ) == 8
	assert median( { 5: 1, 8: 4 } ) == 8
	assert median( { 5: 2, 8: 3 } ) == 8
	assert median( { 5: 3, 8: 2 } ) == 5
	assert median( { 5: 4, 8: 1 } ) == 5
	assert median( { 5: 5, 8: 0 } ) == 5

	# six values
	assert median( { 5: 0, 8: 6 } ) == 8
	assert median( { 5: 1, 8: 5 } ) == 8
	assert median( { 5: 2, 8: 4 } ) == 8
	assert median( { 5: 3, 8: 3 } ) == 5 # see doc
	assert median( { 5: 4, 8: 2 } ) == 5
	assert median( { 5: 5, 8: 1 } ) == 5
	assert median( { 5: 6, 8: 0 } ) == 5
