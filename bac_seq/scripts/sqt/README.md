sqt - SeQuencing Tools
======================

*sqt* is a collection of useful command-line tools for working with high-throughput sequencing data.

Each *sqt* command is a seperate binary program or script with the prefix `sqt-`.

This architecture allows each command to be implemented in any programming language. Simple, one-off scripts can be written in a high-level scripting language such as Python or Perl and later, when performance turns out to be critical, be converted to fast binaries written in a compiled language such as C or C++.

Many *sqt* subcommands are currently implemented in Python. For them, a Python package is available with functions for reading and writing FASTA/FASTQ files, computing alignments, quality trimming, etc.

We welcome submission of new tools! Since on a technical level, there is almost no defined API, except that subcommands simply need to be callable binaries, all subcommands must observe some guidelines in order to offer a consistent interface. For example, the exit code must be zero on success and each tool must offer a `--help` command-line parameter.


Project homepage
----------------

<https://bitbucket.org/marcelm/sqt>


License
-------

(This is the so-called MIT or X11 license.)

* Copyright (c) 2009-2014 Marcel Martin <marcel.martin@scilifelab.se>
* Copyright (c) 2010,2011 Tobias Marschall <tobias.marschall@tu-dortmund.de>
* Copyright (c) 2011 Sven Rahmann <sven.rahmann@tu-dortmund.de>
* Copyright (c) 2012-2013 Johannes Köster <johannes.koester@tu-dortmund.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


Dependencies
------------

- Python
- [cutadapt](https://github.com/marcelm/cutadapt)
- Pysam


Preliminary Guidelines
----------------------

All tools (subcommands) must:

* accept a `--help` parameter that displays help
* give useful exit codes: 0 on success, nonzero on error

The sqt python package is planned to:

* be compatible with Python 2.6 and higher, including Python 3
* be backwards-compatible with older versions of sqt.
  Backwards-compatibility is measured by running the tests:
  If they pass without changes to the tests, everything is ok.


List of Tools
-------------

`sqt-coverage` -- Compute per-reference statistics such as coverage and GC content

`sqt-fastqmod` -- FASTQ modifications: shorten, subset, reverse complement, quality trimming.

`sqt-fastastats` -- Compute N50, min/max length, GC content etc. of a FASTA file

`sqt-qualityguess` -- Guess quality encoding of one or more FASTA files.

`sqt-globalalign` -- Compute a global or semiglobal alignment of two strings.

`sqt-chars` -- Count length of the first word given on the command line.

`sqt-sam-cscq` -- Add the CS and CQ tags to a SAM file with colorspace reads.

`sqt-fastamutate` -- Add substitutions and indels to sequences in a FASTA file.

`sqt-fastaextract` -- Efficiently extract one or more regions from an indexed FASTA file.

`sqt-translate` -- Replace characters in FASTA files (like the 'tr' command).

`sqt-sam-fixn` -- Replace all non-ACGT characters within reads in a SAM file.

`sqt-sam-insertsize` -- Mean and standard deviation of paired-end insert sizes.

`sqt-sam-set-op` -- Set operations (union, intersection, ...) on SAM/BAM files.

`sqt-bam-eof` -- Check for the End-Of-File marker in compressed BAM files.

`sqt-checkfastqpe` -- Check whether two FASTQ files contain correctly paired paired-end data.
