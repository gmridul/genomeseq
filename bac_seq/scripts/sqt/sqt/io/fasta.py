"""
FASTA I/O
"""
from __future__ import print_function
from collections import namedtuple, OrderedDict, Counter
from itertools import islice
import mmap
import sys
from os.path import splitext
from .xopen import xopen
from ..compat import force_str

__author__ = "Marcel Martin"

# use this constant with bytes.translate to convert encoded quality values from
# ascii(phred_quality + 64) to ascii(phred_quality + 33)
TRANSLATE_64_TO_33 = ''.join(chr(max(c-64+33,0)) for c in range(256))


def _shorten(s, n=20):
	"""Shorten string s to at most n characters, appending "..." if necessary."""
	if s is None:
		return None
	if len(s) > n:
		s = s[:n-3] + '...'
	return s


def _quality_to_ascii(qualities, base=33):
	"""
	Convert a list containing qualities given as integer to a string of
	ASCII-encoded qualities.

	base -- ASCII code of quality zero (sensible values are 33 and 64).

	>>> _quality_to_ascii([17, 4, 29, 18])
	'2%>3'
	"""
	qualities = ''.join(chr(q+base) for q in qualities)
	return qualities


class UnknownFileType(Exception):
	"""
	Raised when SequenceReader could not autodetect the file type.
	"""
	pass


class FastaWriter:
	"""
	Write FASTA-formatted sequences to a file-like object.
	"""
	def __init__(self, file, line_length=None):
		"""
		If line_length is not None, the lines will
		be wrapped after line_length characters.
		"""
		self.line_length = line_length
		if isinstance(file, str):
			file = xopen(file, "w")
		self._file = file

	def write(self, name_or_seq, sequence=None):
		"""Write an entry to the the FASTA file.

		If only one parameter (name_or_seq) is given, it must have
		attributes .name and .sequence, which are then used.
		Otherwise, the first parameter must be the name and the second
		the sequence.

		The effect is that you can write this:
		fr.write("name", "ACCAT")
		or
		fr.write(Sequence("name", "ACCAT"))
		"""
		if sequence is None:
			name = name_or_seq.name
			sequence = name_or_seq.sequence
		else:
			name = name_or_seq
		if self.line_length is not None:
			print('>{}'.format(name), file=self._file)
			for i in range(0, len(sequence), self.line_length):
				print(sequence[i:i+self.line_length], file=self._file)
		else:
			print('>{}'.format(name), sequence, file=self._file, sep='\n')

	def close(self):
		self._file.close()

	def __enter__(self):
		if self._file.closed:
			raise ValueError("I/O operation on closed file")
		return self

	def __exit__(self, *args):
		self.close()


class FastqWriter:
	"""
	Write sequences with qualities in FASTQ format.

	FASTQ files are formatted like this:
	@read name
	SEQUENCE
	+
	QUALITIS
	"""
	def __init__(self, file, twoheaders=False):
		"""
		If twoheaders is set, then the read name will be repeated after
		the plus sign (which is redundant and therefore not
		recommended).
		"""
		self.twoheaders = twoheaders
		if isinstance(file, str):
			file = xopen(file, "w")
		self._file = file

	def write(self, name_or_seq, sequence=None, qualities=None):
		"""Write an entry to the the FASTQ file.

		If only one parameter (name_or_seq) is given, it must have
		attributes .name, .sequence and .qualities, which are then
		used. Otherwise, all three parameters must be given and
		name_or_seq must be the name of the sequence.

		The effect is that you can write this:
		fq.write("name", "ACCAT", "#!!&B")
		or
		fq.write(Sequence("name", "ACCAT", "#!!&B"))
		"""
		if sequence is None:
			name = name_or_seq.name
			sequence = name_or_seq.sequence
			qualities = name_or_seq.qualities
		else:
			name = name_or_seq
		if self.twoheaders:
			two = name
		else:
			two = ''
		print("@{0}\n{1}\n+{2}\n{3}".format(name, force_str(sequence), two, force_str(qualities)), file=self._file)

	def close(self):
		self._file.close()

	def __enter__(self):
		if self._file.closed:
			raise ValueError("I/O operation on closed file")
		return self

	def __exit__(self, *args):
		self.close()


IndexEntry = namedtuple("IndexEntry", "name length offset nucleotides_per_line bytes_per_line")


def _read_fasta_index(path):
	"""
	Return a dictionary that maps sequence names to IndexEntry tuples.

	The sequence name is the string before the first space character in the FASTA comment header.

	The offsets are converted to Python 0-based coordinates!
	"""
	index = OrderedDict()
	with open(path) as f:
		for line in f:
			name, length, offset, nucleotides_per_line, bytes_per_line = line.split('\t')
			indexed_name = name.split(' ', 1)[0]
			length = int(length)
			offset = int(offset)
			nucleotides_per_line = int(nucleotides_per_line)
			bytes_per_line = int(bytes_per_line)
			entry = IndexEntry(name, length, offset, nucleotides_per_line, bytes_per_line)
			index[indexed_name] = entry
	return index


class IndexedSequence:
	"""A single sequence in an indexed FASTA file"""

	def __init__(self, mapped, indexentry):
		self.mapped = mapped
		self.indexentry = indexentry

	def __getitem__(self, key):
		if type(key) is int:
			key = slice(key, key+1)
		assert key.step is None

		return self.read(key.start, key.stop)

	def __len__(self):
		"""Return length of sequence"""
		return self.indexentry.length

	def read(self, start=0, stop=None):
		"""
		Read and return a substring of a specific entry of the FASTA file
		"""
		indexinfo = self.indexentry
		if stop is None:
			stop = indexinfo.length
		if start is None:
			start = 0
		start = max(0, start)
		stop = min(indexinfo.length, stop)
		if stop <= start:
			return b''

		def nucleotide_to_byte_offset(i):
			return indexinfo.offset + i // indexinfo.nucleotides_per_line * indexinfo.bytes_per_line + i % indexinfo.nucleotides_per_line

		byte_start = nucleotide_to_byte_offset(start)
		byte_stop = nucleotide_to_byte_offset(stop - 1) + 1

		# these indices are relative to this sequence
		start_line = start // indexinfo.nucleotides_per_line
		stop_line = (stop - 1) // indexinfo.nucleotides_per_line + 1

		assert start_line < stop_line

		if start_line + 1 == stop_line:
			# requested substring is within one line
			return self.mapped[byte_start:byte_stop]

		# otherwise, collect lines or line fragments
		collected = []
		# offset of end of start_line
		byte_line_stop = indexinfo.offset + start_line * indexinfo.bytes_per_line + indexinfo.nucleotides_per_line
		assert byte_start < byte_line_stop
		collected.append(self.mapped[byte_start:byte_line_stop])

		# all lines between the start_- and stop_line
		for line in range(start_line + 1, stop_line - 1):
			byte_line_start = indexinfo.offset + line * indexinfo.bytes_per_line
			assert byte_line_start < byte_line_start + indexinfo.nucleotides_per_line
			collected.append(self.mapped[byte_line_start:byte_line_start+indexinfo.nucleotides_per_line])

		# last line
		byte_line_start = indexinfo.offset + (stop_line - 1) * indexinfo.bytes_per_line
		assert byte_line_start < byte_stop
		collected.append(self.mapped[byte_line_start:byte_stop])

		sequence = b''.join(collected)
		assert len(sequence) == stop - start
		return sequence


class IndexedFasta:
	"""Efficient access to FASTA files that have been indexed by 'samtools faidx'"""

	def __init__(self, path):
		self.index = _read_fasta_index(path + '.fai')
		self.fp = open(path)
		self.mapped = mmap.mmap(self.fp.fileno(), 0, flags=mmap.MAP_PRIVATE) # 0: entire file

	def length(self, name):
		"""Return length of sequence 'name'. TODO make sth. like len(indexedfasta[name]) work"""
		return self.index[name].length

	def get(self, name):
		"""
		Retrieve an entry of the FASTA file as an IndexedSequence object.
		This is a light-weight operation: The actual work happens when you
		access the returned object.
		"""
		return IndexedSequence(self.mapped, self.index[name])

	def __len__(self):
		return len(self.index)

	def __contains__(self, name):
		return name in self.index

	def close(self):
		self.mapped.close()
		self.fp.close()

	def __enter__(self):
		if self.fp.closed:
			raise ValueError("I/O operation on closed file")
		return self

	def __exit__(self, *args):
		self.close()



class Sequence:
	"""qualities is a string and it contains the qualities encoded as ascii(qual+33)."""

	def __init__(self, name, sequence, qualities=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities

	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(self.name, self.sequence[key], self.qualities[key] if self.qualities is not None else None)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = '\', qualities="{0}"'.format(_shorten(self.qualities))
		return 'Sequence(name="{0}", sequence="{1}"{2})'.format(_shorten(self.name), _shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __eq__(self, other):
		return self.name == other.name and \
			self.sequence == other.sequence and \
			self.qualities == other.qualities

	def __ne__(self, other):
		return not self.__eq__(other)


class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""
	pass


class FileWithPrependedLine(object):
	"""
	A file-like object that allows to "prepend" a single
	line to an already opened file. That is, further
	reads on the file will return the provided line and
	only then the actual content. This is needed to solve
	the problem of autodetecting input from a stream:
	As soon as the first line has been read, we know
	the file type, but also that line is "gone" and
	unavailable for further processing.
	"""
	def __init__(self, file, line):
		"""
		file is an already opened file-like object.
		line is a single string (newline will be appended if not included)
		"""
		if not line.endswith('\n'):
			line += '\n'
		self.first_line = line
		self.file = file

	def __iter__(self):
		yield self.first_line
		for line in self.file:
			yield line



def SequenceReader(file, colorspace=False, fileformat=None):
	"""
	Reader for FASTA and FASTQ files that autodetects the file format.
	Returns either an instance of FastaReader or of FastqReader,
	depending on file type.

	The autodetection can be skipped by setting fileformat to the string
	'fasta' or 'fastq'

	file is a filename or a file-like object.
	If file is a filename, then .gz files are supported.
	If the file name is available, the file type is detected
	by looking at the file name.
	If the file name is not available (for example, reading
	from standard input), then the file is read and the file
	type determined from the content.
	"""
	if fileformat is not None:
		fileformat = fileformat.lower()
		if fileformat == 'fasta':
			return FastaReader(file)
		elif fileformat == 'fastq':
			return FastqReader(file)
		else:
			raise UnknownFileType("File format {0} is unknown (expected 'fasta' or 'fastq').".format(fileformat))

	name = None
	if file == "-":
		file = sys.stdin
	elif isinstance(file, str):
		name = file
	elif hasattr(file, "name"):
		name = file.name
	if name is not None:
		if name.endswith('.gz'):
			name = name[:-3]
		name, ext = splitext(name)
		ext = ext.lower()
		if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
			return FastaReader(file)
		elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
			return FastqReader(file, colorspace)
		else:
			raise UnknownFileType("Could not determine whether this is FASTA or FASTQ: file name extension {0} not recognized".format(ext))

	# No name available.
	# Assume that 'file' is an open file
	# and autodetect its type by reading from it.
	for line in file:
		if line.startswith('#'):
			# Skip comment lines (needed for csfasta)
			continue
		if line.startswith('>'):
			return FastaReader(FileWithPrependedLine(file, line))
		if line.startswith('@'):
			return FastqReader(FileWithPrependedLine(file, line), colorspace)
	raise UnknownFileType("File is neither FASTQ nor FASTA.")



class FastaReader(object):
	"""
	Reader for FASTA files.
	"""
	def __init__(self, file, wholefile=False, keep_linebreaks=False, binary=False):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.
		If wholefile is True, then it is ok to read the entire file
		into memory. This is faster when there are many newlines in
		the file, but may obviously need a lot of memory.
		keep_linebreaks -- whether to keep the newline characters in the sequence
		binary -- In Python 3, when this is set, the returned Sequence objects
		 will contain bytes objects for the sequence. The names will still be strings.
		 If file was given as a file-like object, it must have been opened in binary mode.
		"""
		if isinstance(file, str):
			file = xopen(file, 'rb' if binary else 'rt')
		self.fp = file
		self.binary = binary
		self.wholefile = wholefile
		self.keep_linebreaks = keep_linebreaks
		assert not (wholefile and keep_linebreaks), "not supported"

	def __iter__(self):
		"""
		Yield Sequence objects.
		The qualities attribute is always None.
		"""
		return self._wholefile_iter() if self.wholefile else self._streaming_iter()

	def _streaming_iter(self):
		"""
		Read next entry from the file (single entry at a time).
		"""
		name = None
		binary = self.binary
		seq = bytearray() if binary else ""
		delim = b"\n" if binary else '\n'
		startchar = ord('>') if binary else '>'
		for line in self.fp:
			# strip() should also take care of DOS line breaks
			line = line.strip()
			if line and line[0] == startchar:
				if name is not None:
					assert self.keep_linebreaks or seq.find(delim) == -1
					if binary:
						name = name.decode('utf-8')
					if binary:
						seq = bytes(seq)
					yield Sequence(name, seq, None)
				name = line[1:]
				seq = bytearray() if binary else ""
			else:
				seq += line
				if self.keep_linebreaks:
					seq += b'\n' if binary else '\n'
		if name is not None:
			assert self.keep_linebreaks or seq.find(delim) == -1
			if binary:
				name = name.decode('utf-8')
			if binary:
				seq = bytes(seq)
			yield Sequence(name, seq, None)

	def _wholefile_iter(self):
		"""
		This reads in the entire file at once, but is faster than the above code when there are lots of newlines.
		The idea comes from the TAMO package (http://fraenkel.mit.edu/TAMO/), module TAMO.seq.Fasta (author is
		David Benjamin Gordon).
		"""
		if self.binary:
			raise NotImplementedError
		wholefile = self.fp.read()
		assert '\r' not in wholefile, "Sorry, currently don't know how to deal with files that contain \\r linebreaks"
		assert len(wholefile) == 0 or wholefile[0] == '>', "FASTA file must start with '>'"
		parts = wholefile.split('\n>')
		# first part has '>' in front
		parts[0] = parts[0][1:]
		for part in parts:
			lines = part.split('\n', 1)
			yield Sequence(name=lines[0], sequence=lines[1].replace('\n', ''), qualities=None)

	def close(self):
		self.fp.close()

	def __enter__(self):
		if self.fp.closed:
			raise ValueError("I/O operation on closed FastaReader")
		return self

	def __exit__(self, *args):
		self.fp.close()


class FastqReader(object):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, colorspace=False, mode='rt'):
		"""
		file is a filename or a file-like object.
		If file is a filename, then the file is opened with xopen().

		colorspace -- Usually (when this is False), there must be n characters in the sequence and
		 n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.

		mode -- For Python 3, set this to 'rb' to get Sequence objects in which
		 both the sequence and the qualities field have type bytes.
		"""
		if mode not in ('rt', 'rb'):
			raise ValueError("mode must be either 'rt' or 'rb'")
		if isinstance(file, str):
			file = xopen(file, mode)
		self.file = file
		self.colorspace = colorspace
		self.twoheaders = False
		self.mode = mode

	def __iter__(self):
		"""
		Yield Sequence objects.
		"""
		if self.mode == 'rt':
			AT = '@'
			PLUS = '+'
			STRIP = '\n\r'
			name_transform = lambda name: name
		else:
			AT = b'@'
			PLUS = b'+'
			STRIP = b'\n\r'
			name_transform = lambda name: name.decode('utf-8')
		lengthdiff = 1 if self.colorspace else 0
		for i, line in enumerate(self.file):
			if i % 4 == 0:
				if not line.startswith(AT):
					raise FormatError("at line {0}, expected a line starting with '@'".format(i+1))
				name = name_transform(line.strip()[1:])
			elif i % 4 == 1:
				sequence = line.strip()
			elif i % 4 == 2:
				line = line.strip()
				if not line.startswith(PLUS):
					raise FormatError("at line {0}, expected a line starting with '+'".format(i+1))
				if len(line) > 1:
					self.twoheaders = True
					if not name_transform(line[1:]) == name:
						raise FormatError(
							"At line {0}: Two sequence descriptions are given in "
							"the FASTQ file, but they don't match "
							"('{1}' != '{2}')".format(i+1, name, line.rstrip()[1:]))
			elif i % 4 == 3:
				qualities = line.rstrip(STRIP)
				if len(qualities) + lengthdiff != len(sequence):
					raise ValueError("Length of quality sequence and length of read do not match (%d+%d!=%d)" % (len(qualities), lengthdiff, len(sequence)))
				yield Sequence(name, sequence, qualities)

	def __enter__(self):
		if self.file is None:
			raise ValueError("I/O operation on closed FastqReader")
		return self

	def __exit__(self, *args):
		self.file.close()


class FastaQualReader(object):
	"""
	Reader for reads that are stored in .(CS)FASTA and .QUAL files.
	"""
	def __init__(self, fastafile, qualfile, colorspace=False):
		"""
		fastafile and qualfile are filenames file-like objects.
		If file is a filename, then .gz files are supported.

		colorspace -- Usually (when this is False), there must be n characters in the sequence and
		n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.
		"""
		self.fastareader = FastaReader(fastafile)
		self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
		self.colorspace = colorspace

	def __iter__(self):
		"""
		Yield Sequence objects.
		"""
		lengthdiff = 1 if self.colorspace else 0
		for fastaread, qualread in zip(self.fastareader, self.qualreader):
			qualities = _quality_to_ascii(list(map(int, qualread.sequence.split())))
			assert fastaread.name == qualread.name
			if len(qualities) + lengthdiff != len(fastaread.sequence):
				raise ValueError("Length of quality sequence and length of read do not match (%d+%d!=%d)" % (
					len(qualities), lengthdiff, len(fastaread.sequence)))
			yield Sequence(fastaread.name, fastaread.sequence, qualities)

	def __enter__(self):
		if self.fastafile is None:
			raise ValueError("I/O operation on closed FastaQualReader")
		return self

	def __exit__(self, *args):
		self.fastareader.close()
		self.qualreader.close()


def byte_frequencies(s):
	return Counter(s)


def guess_quality_base(path, limit=10000):
	"""Guess quality encoding (Phred33/64) of a FASTQ file.

	Return a tuple (freqs, guess) where freqs are the character frequencies
	encountered in the first limit records and guess is one of 33, 64, None,
	indicating which quality encoding this file probably has (None indicates
	unknown).
	"""
	freqs = Counter()
	with FastqReader(path, mode='rb') as fqr:
		for record in islice(fqr, 0, limit):
			freqs += byte_frequencies(record.qualities)
	if min(freqs) - 64 < -10:
		guess = 33
	elif max(freqs) - 33 > 60:
		guess = 64
	else:
		guess = None
	return freqs, guess


try:
	from sqt._helpers import byte_frequencies
except:
	pass
