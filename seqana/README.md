# Tunneling and sequence analysis
This directory contains code to demonstrate the use of tunneling within the field of sequence analysis.

## Requirements
- GNU make
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [pdflatex](https://linux.die.net/man/1/pdflatex)

## What is contained
- A bundle of programs and a benchmark for the construction of a tunneled FM-index using de Bruijn graph edge minimization.
- A bundle of programs and a benchmark for the construction of a trie, as well as a performant tool to use (compressed) tries
  for multi-pattern search in the Aho-Corasick algorithm.

## Program compilation
To compile all programs, just call `make`.
The software uses the Succinct Data Structure Library [sdsl-lite](https://github.com/simongog/sdsl-lite) by Simon Gog.

After installing sdsl lite, you can either
- set up an environment variable called `SDSLLITE` that specifies the path to the library
  using the command `export SDSLLITE="path/to/sdsllite/"`.
- modify the file `Make.helper` in the root directory such that the uppermost path points to the Make.helper file
  of your sdsl-lite installation.

## De Bruijn graph edge minimization and tunneled FM-index construction

### What is contained
- A program `tfm_index_construct.x` to construct a tunneled FM-index from a file.
  The file must not contain nullbytes, further information can be found by executing the program without parameter.
- A program `tfm_index_invert.x` which can be used to recover the original string from which the FM-index was built.
  Further information can be found by executing the program without parameter.
- A program `dbg_edgespectrum.x` which lists the amount of edges of an edge-reduced de Bruijn graph for the given file
  and a range of different de Bruijn graph orders. The file must not contain nullbytes.
  Further information can be found by executing the program without parameter.

### Benchmark
To execute the benchmark, change to the `benchmark-dbg` directory and call `make`. The benchmark can be configured using the file
`benchmark.config`. The benchmark should produce a pdf file containing the results, similar to [this one](benchmark-dbg/dbgbenchmark.pdf).

Depending on the size of the test files, it is useful to have a platform with at least 16 GB of memory.

## Trie construction and multi-pattern search using the Aho-Corasick algorithm

### What is contained
- A program `create_trie_input.x` to prepare a file of newline-separated strings for the construction of a trie.
  The program removes nullbytes from the strings and ensures that no string is a proper substring of another string,
  to ensure correctness of the multi-pattern search with the Aho-Corasick algorithm.
  Further information can be found by executing the program without parameter.
- A program `trie_construct.x` which can be used to construct a trie from the given input file.
  The input file must not contain nullbytes. The lines of the file are used as inputs for the trie.
  No line should be a proper substring of another file to ensure the Aho-Corasick algorithm works correctly.
  Further information can be found by executing the program without parameter.
- A program `mp_search.x` which searches for all strings in the given trie within the standard input.
  The trie should be created using the `create_trie_input.x` and the `trie_construct.x` programs to ensure correctness,
  because the implemented Aho-Corasick algorithm is a simplified version.

### Benchmark
To execute the benchmark, change to the `benchmark-trie` directory and call `make`. The benchmark can be configured using the file
`benchmark.config`. The benchmark should produce a pdf file containing the results, similar to [this one](benchmark-trie/triebenchmark.pdf).
Note that the benchmark filters lines with less than 10 characters from the output.

Depending on the size of the test files, it is useful to have a platform with at least 16 GB of memory.
