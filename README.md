# sequence-analysis

Programs for analyzing biological sequences such as DNA, RNA, and proteins.

This project currently implements Needleman-Wunsch and
Waterman-Smith-Beyer algorithms for pairwise sequence alignment.

Reports are generated using Sphinx Doc.


# Prerequisites

- Python 3
- Linux (might work elsewhere too)

Building PDFs also requires a full LaTeX installation.


# Usage

To install (this will create and install everything into a Python
virtualenv):

```shell
$ ./install.sh
```

To run the tests and code quality checks:
```shell
$ ./run-tests.sh
```

For algorithms, please see the doc:
https://github.com/ajsmith/sequence-analysis/blob/main/hw1/source/index.rst

For tests, see the "tests/" directory.
