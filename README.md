# miRador

[![Build Status](https://travis-ci.org/PrashntS/pancake.svg?branch=master)](https://travis-ci.org/PrashntS/pancake)

miRador finds potential target sites for the microRNAs in genomic sequences. It is built on miRanda, an algorithm for detection and ranking of the targets of microRNA.

## Usage

## Installation

## Description

The miRanda algorithm works in two phases. In phase one, the potential target sites are reported based on query microRNA and reference (CDNA) sequence. These targets are scored and the high scoring alignments are then used in second phase, where the folding routines of RNAlib library are utilised to calculate the minimum free energy of the resulting combinations.

## Further Information

- miRanda algorithm on [microrna.org](http://www.microrna.org/microrna/getDownloads.do) | [user manual](http://cbio.mskcc.org/microrna_data/manual.html)
- RNAlib at [Vienna RNA](http://www.tbi.univie.ac.at/RNA/)

### Citing in publications

Please cite the original miRanda library, and Vienna RNA library. The citations can be obtained from the links above.
