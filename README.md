# Outrigger

[![Build Status](https://travis-ci.org/YeoLab/outrigger.svg?branch=master)](https://travis-ci.org/YeoLab/outrigger)[![](https://img.shields.io/pypi/v/outrigger.svg)](https://pypi.python.org/pypi/outrigger)[![Coverage Status](https://coveralls.io/repos/YeoLab/outrigger/badge.svg?branch=master&service=github)](https://coveralls.io/github/YeoLab/outrigger?branch=master)

Outrigger is a program to calculate splicing scores of RNA-Seq data based on junction reads and a *de novo*, custom annotation created with a graph database.

 1. Read all your `SJ.out.tab` files from the STAR aligner into a single, compiled file
 2. Building a *de novo* splicing annotation index [of skipped exon (SE) and mutually exclusive exon (MXE) events] that is custom to the junctions observed in your data
 3. (optional) Consolidate events that share junctions but have alternative start and ends for the flanking exons. The criteria to consolidate these is based on which pair of isoforms contains one isoform that is annotated as the "principal" isoform by APPRIS (in the GENCODE annotation).
 4. Calculate "percent spliced-in" (Psi/Ψ) scores for all your samples given the events calculated in (2)

 The program `outrigger index` takes care of (1-3) for you. `outrigger psi` does step (4). If you're using an existing index with new junction reads, then `outrigger psi` will also do (1).

* Free software: BSD license

## Installation

### Bleeding edge code from Github (here)

For advanced users, if you have [git](https://git-scm.com/) and [Anaconda Python](https://www.continuum.io/downloads) installed, you can:

1. Clone this repository
2. Change into that directory
3. Create an environment with the necessary packages from Anaconda
4. Activate the environment
5. Install remaining packages from PyPI ([`graphlite`](https://github.com/eugene-eeo/graphlite) is only available on PyPI, not as a `conda` package)
6. Install this package

These steps are shown in code below.

```
git clone git@github.com:YeoLab/outrigger
cd outrigger
conda create --name outrigger --yes --file conda_requirements.txt --channel bioconda
source activate outrigger
pip install -r requirements.txt
pip install .
```


## Features

### `index`: Build a *de novo* splicing annotation index of events custom to *your* data

Help features

```
usage: outrigger.py index [-h] [-i INDEX] -j [SJ_OUT_TAB [SJ_OUT_TAB ...]]
                          [-m MIN_READS] (-g GTF | -d GFFUTILS_DATABASE)

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        Name of the folder where you saved the output from
                        "outrigger index" (default is ./outrigger_index, which
                        is relative to the directory where you called the
                        program)". You will need this file for the next step,
                        "outrigger psi"
  -j [SJ_OUT_TAB [SJ_OUT_TAB ...]], --sj-out-tab [SJ_OUT_TAB [SJ_OUT_TAB ...]]
                        SJ.out.tab files from STAR aligner output
  -m MIN_READS, --min-reads MIN_READS
                        Minimum number of reads per junction for that junction
                        to count in creating the index of splicing events
  -g GTF, --gtf GTF     Name of the gtf file you want to use. If a gffutils
                        feature database doesn't already exist at this
                        location plus '.db' (e.g. if your gtf is
                        gencode.v19.annotation.gtf, then the database is
                        inferred to be gencode.v19.annotation.gtf.db), then a
                        database will be auto-created. Not required if you
                        provide a pre-built database with '--gffutils-
                        database'
  -d GFFUTILS_DATABASE, --gffutils-database GFFUTILS_DATABASE
                        Name of the gffutils database file you want to use.
                        The exon IDs defined here will be used in the function
                        when creating splicing event names. Not required if
                        you provide a gtf file with '--gtf'
```

#### Example command

```
outrigger index --sj-out-tab *SJ.out.tab --gtf gencode.v19.annotation.gtf
```

This is equivalent to the below command, which specifies all the other arguments with the default values:

```
outrigger index --index ./outrigger_index --sj-out-tab *SJ.out.tab --min-reads 10 --gtf gencode.v19.annotation.gtf
```

#### Outputs

The above commands will create a folder called `outrigger_index` in the folder you ran the command from, with the following structure

```
outrigger_index/
    events/
        se.csv
        mxe.csv
    sj.csv
```


### `psi`: Calculate percent spliced-in (Psi/Ψ) scores for your data from the splicing events you created

```
$ outrigger psi --help
usage: outrigger.py psi [-h] -i INDEX
                        [-c SPLICE_JUNCTION_CSV | -j [SJ_OUT_TAB [SJ_OUT_TAB ...]]]
                        [-m MIN_READS] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        Name of the folder where you saved the output from
                        "outrigger index" (default is ./outrigger_index, which
                        is relative to the directory where you called the
                        program)"
  -c SPLICE_JUNCTION_CSV, --splice-junction-csv SPLICE_JUNCTION_CSV
                        Name of the splice junction files to calculate psi
                        scores on. If not provided, the compiled 'sj.csv' file
                        with all the samples from the SJ.out.tab files that
                        were used during 'outrigger index' will be used. Not
                        required if you specify SJ.out.tab file with '--sj-
                        out-tab'
  -j [SJ_OUT_TAB [SJ_OUT_TAB ...]], --sj-out-tab [SJ_OUT_TAB [SJ_OUT_TAB ...]]
                        SJ.out.tab files from STAR aligner output. Not
                        required if you specify
  -m MIN_READS, --min-reads MIN_READS
                        Minimum number of reads per junction for calculating
                        Psi
  --debug               If given, print debugging logging information to
                        standard out
```


```
outrigger psi
```

The above command is equivalent to specifying all the arguments with their default values:

```
outrigger psi --index ./outrigger_index --min-reads 10
```


## Outputs

```
outrigger/
├── index
│   └── mxe
│       ├── junctions.csv
│       └── metadata.csv
│   └── se
│       ├── junctions.csv
│       └── metadata.csv
├── psi.csv
└── junction_reads.csv
```
