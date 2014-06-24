"""
Tests for the sqt.intervaltree module
"""
from nose.tools import raises
import sys
from sqt.intervaltree import IntervalTree

__author__ = "Johannes Köster"

def test_intervaltree():
	tree = IntervalTree()
	tree.insert(1,10)
	tree.insert(5,20)
	tree.insert(8,20)
	print(tree, file=sys.stderr)
	print(list(tree.find(10,11)))
	assert len(list(tree.find(10,11))) == 3
	assert len(list(tree.find(1,4))) == 1
	assert len(list(tree.find(30,40))) == 0
