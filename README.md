# Outrigger

![Outrigger logo](https://raw.githubusercontent.com/YeoLab/outrigger/master/logo/logo_v1.png)

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

The "help" output of the two programs tries to be explicit about what is required
to run `outrigger`. Below is the output of when you use the command,
`outrigger index --help`

```
$ outrigger index --help
usage: outrigger index [-h] [-o OUTPUT]
                       (-j [SJ_OUT_TAB [SJ_OUT_TAB ...]] | -c JUNCTION_READ_CSV)
                       [-m MIN_READS] [--use-multimapping]
                       (-g GTF_FILENAME | -d GFFUTILS_DB) [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Name of the folder where you saved the output from
                        "outrigger index" (default is ./outrigger_output,
                        which is relative to the directory where you called
                        the program)". You will need this file for the next
                        step, "outrigger psi" (default="./outrigger_output")
  -j [SJ_OUT_TAB [SJ_OUT_TAB ...]], --sj-out-tab [SJ_OUT_TAB [SJ_OUT_TAB ...]]
                        SJ.out.tab files from STAR aligner output
  -c JUNCTION_READ_CSV, --junction-read-csv JUNCTION_READ_CSV
                        Name of the splice junction files to calculate psi
                        scores on. If not provided, the compiled
                        './outrigger_output/junction_reads/reads.csv' file
                        with all the samples from the SJ.out.tab files that
                        were used during 'outrigger index' will be used. Not
                        required if you specify SJ.out.tab file with '--sj-
                        out-tab'
  -m MIN_READS, --min-reads MIN_READS
                        Minimum number of reads per junction for that junction
                        to count in creating the index of splicing events
                        (default=10)
  --use-multimapping    Applies to STAR SJ.out.tab files only. If this flag is
                        used, then include reads that mapped to multiple
                        locations in the genome, not uniquely to a locus, in
                        the read count for a junction. By default, this is
                        off, and only uniquely mapped reads are used.
  -g GTF_FILENAME, --gtf-filename GTF_FILENAME
                        Name of the gtf file you want to use. If a gffutils
                        feature database doesn't already exist at this
                        location plus '.db' (e.g. if your gtf is
                        gencode.v19.annotation.gtf, then the database is
                        inferred to be gencode.v19.annotation.gtf.db), then a
                        database will be auto-created. Not required if you
                        provide a pre-built database with '--gffutils-db'
  -d GFFUTILS_DB, --gffutils-db GFFUTILS_DB
                        Name of the gffutils database file you want to use.
                        The exon IDs defined here will be used in the function
                        when creating splicing event names. Not required if
                        you provide a gtf file with '--gtf-filename'
  --debug               If given, print debugging logging information to
                        standard out (Warning: LOTS of output. Not recommended
                        unless you think something is going wrong)
```

#### Example command

Included in this repository is a subset of the 1809 cells from
["Adult mouse cortical cell taxonomy revealed by single cell transcriptomics."
by Tasic et al, Nature Neuroscience (2016)](http://www.ncbi.nlm.nih.gov/pubmed/26727548).
There splice junction output files from the [STAR aligner](https://github.com/alexdobin/STAR) from the 43 "`CAV_LP_Ipsi_tdTpos`" cells,
plus a subset of the [GENCODE M10](http://www.gencodegenes.org/mouse_releases/10.html)
(Version M10 (January 2016 freeze, GRCm38) - Ensembl 85) mouse annotation.

To run this program with the included example data, from the `outrigger` directory
where you cloned `outrigger` (this is important because the locations of the files
is relative to that directory), run this command:

```
outrigger index \
    --sj-out-tab test_data/tasic2016/unprocessed/sj_out_tab/* \
    --gtf test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf
```

*Note: the backslashes (`\`, like a tree that's falling backwards relative to
right-to-left reading) to tell the computer that you're not
done tellingn it what to do, so the line continues, and to help the code be a
little more human-readable. The above code is read by the computer exactly the
same as the one-liner below, but is easier for us humans to read.*

```
outrigger index --sj-out-tab test_data/tasic2016/unprocessed/sj_out_tab/* --gtf test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf
```

This is equivalent to the below command, which specifies all the other arguments
with the default values.

```
outrigger index \
    --sj-out-tab test_data/tasic2016/unprocessed/sj_out_tab/* \
    --gtf test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf \
    --output ./outrigger_output --min-reads 10
```

The output of this command is:


```
$ outrigger index --sj-out-tab example_data/tasic2016/unprocessed/sj_out_tab/* --gtf example_data/tasic2016/unprocessed/gtf/snap25_myl6.gtf
2016-08-12 11:24:03	Reading SJ.out.files and creating a big splice junction table of reads spanning exon-exon junctions...
2016-08-12 11:24:03	Writing ./outrigger_output/junction_reads/reads.csv ...

2016-08-12 11:24:03		Done.
2016-08-12 11:24:03	Creating splice junction metadata of merely where junctions start and stop
2016-08-12 11:24:03		Done.
2016-08-12 11:24:03	Getting junction-direction-exon triples for graph database ...
2016-08-12 11:24:03	Starting annotation of all junctions with known neighboring exons ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Writing junction-exon-direction triples to ./outrigger_output/index/junction_exon_direction_triples.csv...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Populating graph database of the junction-direction-exon triples ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Finding all skipped exon (SE) events ...
2016-08-12 11:24:04	Trying out 25 exons ...
2016-08-12 11:24:04		1/25 exons tested (4.0%)
2016-08-12 11:24:04		2/25 exons tested (8.0%)
2016-08-12 11:24:04		3/25 exons tested (12.0%)
2016-08-12 11:24:04		4/25 exons tested (16.0%)
2016-08-12 11:24:04		5/25 exons tested (20.0%)
2016-08-12 11:24:04		6/25 exons tested (24.0%)
2016-08-12 11:24:04		7/25 exons tested (28.0%)
2016-08-12 11:24:04		8/25 exons tested (32.0%)
2016-08-12 11:24:04		9/25 exons tested (36.0%)
2016-08-12 11:24:04		10/25 exons tested (40.0%)
2016-08-12 11:24:04		11/25 exons tested (44.0%)
2016-08-12 11:24:04		12/25 exons tested (48.0%)
2016-08-12 11:24:04		13/25 exons tested (52.0%)
2016-08-12 11:24:04		14/25 exons tested (56.0%)
2016-08-12 11:24:04		15/25 exons tested (60.0%)
2016-08-12 11:24:04		16/25 exons tested (64.0%)
2016-08-12 11:24:04		17/25 exons tested (68.0%)
2016-08-12 11:24:04		18/25 exons tested (72.0%)
2016-08-12 11:24:04		19/25 exons tested (76.0%)
2016-08-12 11:24:04		20/25 exons tested (80.0%)
2016-08-12 11:24:04		21/25 exons tested (84.0%)
2016-08-12 11:24:04		22/25 exons tested (88.0%)
2016-08-12 11:24:04		23/25 exons tested (92.0%)
2016-08-12 11:24:04		24/25 exons tested (96.0%)
2016-08-12 11:24:04		25/25 exons tested (100.0%)
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Writing 1 SE events to ./outrigger_output/index/se/events.csv ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Making metadata file of SE events, annotating them with GTF attributes ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Getting exon and intron lengths of alternative events ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Combining lengths and attributes into one big dataframe ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Writing SE metadata to ./outrigger_output/index/se/metadata.csv ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Finding all mutually exclusive exon (MXE) events ...
2016-08-12 11:24:04	Trying out 25 exons ...
2016-08-12 11:24:04		1/25 exons tested (4.0%)
2016-08-12 11:24:04		2/25 exons tested (8.0%)
2016-08-12 11:24:04		3/25 exons tested (12.0%)
2016-08-12 11:24:04		4/25 exons tested (16.0%)
2016-08-12 11:24:04		5/25 exons tested (20.0%)
2016-08-12 11:24:04		6/25 exons tested (24.0%)
2016-08-12 11:24:04		7/25 exons tested (28.0%)
2016-08-12 11:24:04		8/25 exons tested (32.0%)
2016-08-12 11:24:04		9/25 exons tested (36.0%)
2016-08-12 11:24:04		10/25 exons tested (40.0%)
2016-08-12 11:24:04		11/25 exons tested (44.0%)
2016-08-12 11:24:04		12/25 exons tested (48.0%)
2016-08-12 11:24:04		13/25 exons tested (52.0%)
2016-08-12 11:24:04		14/25 exons tested (56.0%)
2016-08-12 11:24:04		15/25 exons tested (60.0%)
2016-08-12 11:24:04		16/25 exons tested (64.0%)
2016-08-12 11:24:04		17/25 exons tested (68.0%)
2016-08-12 11:24:04		18/25 exons tested (72.0%)
2016-08-12 11:24:04		19/25 exons tested (76.0%)
2016-08-12 11:24:04		20/25 exons tested (80.0%)
2016-08-12 11:24:04		21/25 exons tested (84.0%)
2016-08-12 11:24:04		22/25 exons tested (88.0%)
2016-08-12 11:24:04		23/25 exons tested (92.0%)
2016-08-12 11:24:04		24/25 exons tested (96.0%)
2016-08-12 11:24:04		25/25 exons tested (100.0%)
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Writing 1 MXE events to ./outrigger_output/index/mxe/events.csv ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Making metadata file of MXE events, annotating them with GTF attributes ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Getting exon and intron lengths of alternative events ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Combining lengths and attributes into one big dataframe ...
2016-08-12 11:24:04		Done.
2016-08-12 11:24:04	Writing MXE metadata to ./outrigger_output/index/mxe/metadata.csv ...
2016-08-12 11:24:04		Done.
```


#### Outputs

The above commands will create a folder called `outrigger_index` in the folder you ran the command from, with the following structure

```
$ tree outrigger_output
outrigger_output
├── index
│   ├── junction_exon_direction_triples.csv
│   ├── mxe
│   │   ├── events.csv
│   │   └── metadata.csv
│   └── se
│       ├── events.csv
│       └── metadata.csv
└── junction_reads
    └── reads.csv

4 directories, 6 files
```


### `psi`: Calculate percent spliced-in (Psi/Ψ) scores for your data from the splicing events you created

```
$ outrigger psi --help
usage: outrigger psi [-h] [-i INDEX]
                     [-c JUNCTION_READ_CSV | -j [SJ_OUT_TAB [SJ_OUT_TAB ...]]]
                     [-m MIN_READS] [--use-multimapping]
                     [--reads-col READS_COL] [--sample-id-col SAMPLE_ID_COL]
                     [--junction-id-col JUNCTION_ID_COL] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        Name of the folder where you saved the output from
                        "outrigger index" (default is ./outrigger_index, which
                        is relative to the directory where you called this
                        program, assuming you have called "outrigger psi" in
                        the same folder as you called "outrigger index")
  -c JUNCTION_READ_CSV, --junction-read-csv JUNCTION_READ_CSV
                        Name of the splice junction files to calculate psi
                        scores on. If not provided, the compiled
                        './outrigger_output/junction_reads/reads.csv' file
                        with all the samples from the SJ.out.tab files that
                        were used during 'outrigger index' will be used. Not
                        required if you specify SJ.out.tab file with '--sj-
                        out-tab'
  -j [SJ_OUT_TAB [SJ_OUT_TAB ...]], --sj-out-tab [SJ_OUT_TAB [SJ_OUT_TAB ...]]
                        SJ.out.tab files from STAR aligner output. Not
                        required if you specify a file with "--junction-read-
                        csv"
  -m MIN_READS, --min-reads MIN_READS
                        Minimum number of reads per junction for calculating
                        Psi (default=10)
  --use-multimapping    Applies to STAR SJ.out.tab files only. If this flag is
                        used, then include reads that mapped to multiple
                        locations in the genome, not uniquely to a locus, in
                        the read count for a junction. By default, this is
                        off, and only uniquely mapped reads are used.
  --reads-col READS_COL
                        Name of column in --splice-junction-csv containing
                        reads to use. (default='reads')
  --sample-id-col SAMPLE_ID_COL
                        Name of column in --splice-junction-csvcontaining
                        sample ids to use. (default='sample_id')
  --junction-id-col JUNCTION_ID_COL
                        Name of column in --splice-junction-csvcontaining the
                        ID of the junction to use. Must match exactly with the
                        junctions in the index.(default='junction_id')
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
