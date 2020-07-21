# Tunneling and data compression
This directory contains code to enhance BWT-based compressors with tunneling.
Moreover, a data compression benchmark is included.

## Requirements
- GNU make
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [pdflatex](https://linux.die.net/man/1/pdflatex)

## What is contained
- A program `bwzip.x` containing different BWT-based compressors which can be enhanced with tunneling.
  Call the program without a parameter to see information on usage.
- A program `tfmzip.x` which can be used to compress tunneled (or normal) FM-indices, see [sequence analysis](../seqana).
  Call the program without a parameter to see information on usage.
- A data compression benchmark. 

## Program compilation
To compile the programs `bwzip.x` and `tfmzip.x`, just call `make`.
The software uses the Succinct Data Structure Library [sdsl-lite](https://github.com/simongog/sdsl-lite) by Simon Gog.

After installing sdsl lite, you can either
- set up an environment variable called `SDSLLITE` that specifies the path to the library
  using the command `export SDSLLITE="path/to/sdsllite/"`.
- modify the file `Make.helper` in the root directory such that the uppermost path points to the Make.helper file
  of your sdsl-lite installation.

## Benchmark
A benchmark is contained under the `benchmark` directory.
The benchmark uses test data from the [test data directory](../testdata/). Beside of BWT-based compressors,
the benchmark includes the compressors [xz](https://tukaani.org/xz/) and [zpaq](http://mattmahoney.net/dc/zpaq.html).

Test data to be used can be configured in the file `benchmark.config`.
The benchmark can be executed as follows. Go into the `benchmark` directory and run the following:
1. call `make install` to install and compile the required programs (you may will need superuser rights for this, so use `sudo make install`).
2. call `make benchmark` to execute the benchmark itself. Depending on the used test files, it might be useful to have 8 GB of memory for files with less than 512 MB size,
   16 GB of memory for files with less than 1 GB size, and up to 32 GB of memory for bigger files.
3. call `make visualize` to visualize the results in form of a pdf.
All steps can be performed at once by just calling `make`.
