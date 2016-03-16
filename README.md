# metaRNA

metaRNA finds potential target sites for the microRNAs in genomic sequences. It is built on miRanda, an algorithm for detection and ranking of the targets of microRNA.

[![Build Status](https://travis-ci.org/PrashntS/metaRNA.svg?branch=master)](https://travis-ci.org/PrashntS/metaRNA)
[![PyPI](https://img.shields.io/pypi/v/metarna.svg)](https://pypi.python.org/pypi/metarna)
[![PyPI](https://img.shields.io/pypi/dw/metarna.svg)](https://pypi.python.org/pypi/metarna)
[![PyPI](https://img.shields.io/pypi/pyversions/metarna.svg)](https://pypi.python.org/pypi/metarna)
[![PyPI](https://img.shields.io/pypi/l/metarna.svg)]()
[![](https://img.shields.io/github/issues-raw/prashnts/metarna.svg)](https://github.com/PrashntS/metaRNA/issues)
[![PyPI](https://img.shields.io/pypi/status/metarna.svg)]()

## Quickstart

```python
from metarna.target_scan import scan, free_energy

gene_sequence = (
    "ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCC"
    "CCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGC"
    "CTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGG"
)

mirna_sequence = "UGGCGAUUUUGGAACUCAAUGGCA"

# Get free Energy value:
delta_g = free_energy(gene_sequence, mirna_sequence)

# Get full targets information:
targets = scan(gene_sequence, mirna_sequence)

```

## Installation

metaRNA supports Python versions 2.7, 3.3, 3.4, and 3.5. It requires the Vienna RNA package which must be installed before installing metaRNA.

After Intalling Vienna RNA package, metaRNA may be installed simply by executing:
```shell
$ pip install metarna
```

metaRNA is currently tested on Mac OSX and Ubuntu, however other Unix based systems should be supported. It isn't tested on Windows yet.

## Running tests

Use of `virtualenv` is assumed and expected.

```shell
$ python setup.py develop # Installs Development Version
$ python -m unittest
```

## Description

The miRanda algorithm works in two phases. In phase one, the potential target sites are reported based on query microRNA and reference (CDNA) sequence. These targets are scored and the high scoring alignments are then used in second phase, where the folding routines of RNAlib library are utilised to calculate the minimum free energy of the resulting combinations.

## Further Information

- miRanda algorithm on [microrna.org](http://www.microrna.org/microrna/getDownloads.do) | [user manual](http://cbio.mskcc.org/microrna_data/manual.html)
- RNAlib at [Vienna RNA](http://www.tbi.univie.ac.at/RNA/)

### Citing in publications

Please cite the original miRanda library, and Vienna RNA library. The citations can be obtained from the links above.
