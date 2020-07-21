# BWT Tunneling
This bundle of programs is part of the thesis

	BWT Tunneling
	Uwe Baier, 2020

The bundle consists of three parts:
1. A suite to handle [test data](testdata/).
2. Programs that show the [application of tunneling in data compression](datacomp/).
3. Programs that show the [application of tunneling in sequence analysis](seqana/).

## General requirements
- GNU make
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [pdflatex](https://linux.die.net/man/1/pdflatex)

## Installation
Every part requires its own programs, so you may read the README of each part.
Apart from that, all parts require the Succinct Data Structure Library [sdsl-lite](https://github.com/simongog/sdsl-lite) by Simon Gog.

After installing sdsl lite, you can either
- set up an environment variable called `SDSLLITE` that specifies the path to the library
  using the command `export SDSLLITE="path/to/sdsllite/"`.
- modify the file `Make.helper` in this directory such that the uppermost path points to the Make.helper file
  of your sdsl-lite installation.
