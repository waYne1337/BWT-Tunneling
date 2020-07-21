# Test Data Suite
This little suite helps you to manage test data.

## Requirements
- GNU make
- [curl](https://curl.haxx.se/) for downloading data
- [gzip](http://www.gzip.org/)
- [bzip2](http://www.bzip.org/)
- [zip](http://infozip.sourceforge.net/)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [pdflatex](https://linux.die.net/man/1/pdflatex)

## Before you start
Genome data must be converted into special formats,
which is handled by the utility program 'twoBitsToFa'.
Before you can download any data, be sure to set the
executable bit of this program, i.e. call
`chmod a+x twoBitsToFa`.

The software uses the Succinct Data Structure Library [sdsl-lite](https://github.com/simongog/sdsl-lite) by Simon Gog.

After installing sdsl lite, you can either
- set up an environment variable called `SDSLLITE` that specifies the path to the library
  using the command `export SDSLLITE="path/to/sdsllite/"`.
- modify the file `Make.helper` in the root directory such that the uppermost path points to the Make.helper file
  of your sdsl-lite installation.

## Test Data
Test data sorted by text corpora can be found in the
file `data.config`. You may include this file to your own
makefile to have references to all files of a certain text corpus.
This suite includes the following text corpora:
- [canterbury corpus](http://www.data-compression.info/Corpora/CanterburyCorpus/index.html)
- [largecanterbury corpus](http://www.data-compression.info/Corpora/CanterburyCorpus/index.html)
- [silesia corpus](http://sun.aei.polsl.pl/~sdeor/index.php?page=silesia)
- [pizzachili corpus](http://pizzachili.dcc.uchile.cl/texts.html)
- [repetitive corpus](http://pizzachili.dcc.uchile.cl/repcorpus.html) (only real texts)
- [genomes](http://hgdownload.soe.ucsc.edu/downloads.html) (only a selection of 3 files)

A file with statistics about the test data can be found at [here](filestat.pdf).

## Usage
- To load all included test files, call `make`.
- To fully load a corpus cpname, use `make cpname`.
- To download only a specific file named fname of corpus cpname, 
  use `make cpname/fname`.
- To clean data of a text corpus cpname, use `make clean-cpname`.
- To clean all data, use `make cleanall`

## Additional Notes
You may wish to include the file data.config to your own Makefile,
so the data paths immediately are present. However, when loading 
data, always use the makefile in this directory and DO NOT include it,
because it uses relative paths (possibly use 
`cd data;make cpname/fname1 cpname/fname2 ...` to load your data 
before using it).

## Utility
A little utility program named `filestat.cpp` is included: it computes
alphabet and file size for a given file. To get detailed file statistics
in form of a PDF file, call `make statistics`.
