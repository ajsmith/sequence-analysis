# sequence-analysis

Programs for analyzing biological sequences such as DNA, RNA, and
proteins.

This project currently implements Needleman-Wunsch and
Waterman-Smith-Beyer algorithms for pairwise sequence alignment.

Reports are generated using Sphinx Doc.


# Prerequisites

- Python 3
- Linux (might work elsewhere too)
- LaTeX (optional if building PDFs is needed)


# Usage

Optionally, everything can be run within the `binftools` container
(separate project) which provides all prerequisites listed above. Just
mount this project within the container and get to work.

```shell
$ cd sequence-analysis
$ docker run -it --rm -v "$(pwd):/sequence-analysis" ajsmith/binftools bash
[me@22922011bda7 ~]$ cd /sequence-analysis
```


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
