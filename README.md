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

To install `outrigger`, we recommend using the
[Anaconda Python Distribution](http://anaconda.org/) and creating an environment.

You'll want to add the [`bioconda`](https://bioconda.github.io/) channel to
make installing [`bedtools`](bedtools.readthedocs.io) and its Python wrapper,
[`pybedtools`](https://daler.github.io/pybedtools/) easy.

```
conda config --add channels r
conda config --add channels bioconda

```

Create an environment called `outrigger-env`

```
conda create -n outrigger-env pandas pybedtools gffutils biopython
```

Now activate that environment and install `outrigger` from PyPI:

```
source activate outrigger-env
pip install outrigger
```

To check that it installed properly, try the command with the help option (`-h`), `outrigger -h`. The output
should look like this:

```
$ outrigger -h
usage: outrigger [-h] {index,validate,psi} ...

Calculate "percent-spliced in" (Psi) scores of alternative splicing on a *de
novo*, custom-built splicing index

positional arguments:
  {index,validate,psi}  Sub-commands
    index               Build an index of splicing events using a graph
                        database on your junction reads and an annotation
    validate            Ensure that the splicing events found all have the
                        correct splice sites
    psi                 Calculate "percent spliced-in" (Psi) values using the
                        splicing event index built with "outrigger index"

optional arguments:
  -h, --help            show this help message and exit
```


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

## Quick start

If you just want to know how to run this on your data with the default
parameters, start here. Let's say you performed your alignment in the folder
called `~/projects/tasic2016/analysis/tasic2016_v1`, and that's where your
`SJ.out.tab` files from the STAR aligner are (they're output into the same
folder as the `.bam` files). First you'll need to change directories to that
folder with `cd`.

```
cd ~/projects/tasic2016/analysis/tasic2016_v1
```

Then you need find all alternative splicing events, which you do by running
`outrigger index` on the splice junction files and the gtf. Here is an example
command:

```
outrigger index --sj-out-tab *SJ.out.tab \
    --gtf /projects/ps-yeolab/genomes/mm10/gencode/m10/gencode.vM10.annotation.gtf
```

Next, you'll want to validate that the splicing events you found follow
biological rules, such as being containing GT/AG (mammalian major spliceosome)
or AT/AC (mammalian minor splicesome) sequences. To do that, you'll need to
provide the genome name (e.g. `mm10`) and the genome sequences. An example
command is below:

```
outrigger validate --genome mm10 \
    --fasta /projects/ps-yeolab/genomes/mm10/GRCm38.primary_assembly.genome.fa
```


Finally, you can calculate percent spliced in (Psi) of your splicing events!
Thankfully this is very easy:

```
outrigger psi
```

It should be noted that ALL of these commands should be performed in the same
directory, so no moving.

### Quick start summary




```
cd ~/projects/tasic2016/analysis/tasic2016_v1
outrigger index --sj-out-tab *SJ.out.tab \
    --gtf /projects/ps-yeolab/genomes/mm10/gencode/m10/gencode.vM10.annotation.gtf
outrigger validate --genome mm10 \
    --fasta /projects/ps-yeolab/genomes/mm10/GRCm38.primary_assembly.genome.fa
outrigger psi
```

This will create a folder called `outrigger_output`, which at the end should
look like this:





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
    --sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* \
    --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf
```

*Note: the backslashes (`\`, like a tree that's falling backwards relative to
right-to-left reading) to tell the computer that you're not
done telling it what to do, so the line continues, and to help the code be a
little more human-readable. The above code is read by the computer exactly the
same as the one-liner below, but is easier for us humans to read.*

```
outrigger index --sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf
```

This is equivalent to the below command, which specifies all the other arguments
with the default values.

```
outrigger index \
    --sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* \
    --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf \
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
│   ├── gtf
│   │   ├── gencode.vM10.annotation.snap25.myl6.gtf
│   │   ├── gencode.vM10.annotation.snap25.myl6.gtf.db
│   │   ├── gencode.vM10.annotation.snap25.myl6.gtf.db.bak
│   │   └── novel_exons.gtf
│   ├── junction_exon_direction_triples.csv
│   ├── mxe
│   │   ├── events.csv
│   │   ├── exon1.bed
│   │   ├── exon2.bed
│   │   ├── exon3.bed
│   │   ├── exon4.bed
│   │   └── metadata.csv
│   └── se
│       ├── events.csv
│       ├── exon1.bed
│       ├── exon2.bed
│       ├── exon3.bed
│       └── metadata.csv
└── junctions
    ├── metadata.csv
    └── reads.csv

5 directories, 18 files
```

### `validate`: Check that the found exons are real

```
outrigger validate -f ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa -g ~/genomes/mm10/mm10.chrom.sizes
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

Now the `outrigger_output` folder has `psi` subfolder, with the MXE and SE
events separate.


```
$ tree outrigger_output
outrigger_output
├── index
│   ├── gtf
│   │   ├── gencode.vM10.annotation.snap25.myl6.gtf
│   │   ├── gencode.vM10.annotation.snap25.myl6.gtf.db
│   │   ├── gencode.vM10.annotation.snap25.myl6.gtf.db.bak
│   │   └── novel_exons.gtf
│   ├── junction_exon_direction_triples.csv
│   ├── mxe
│   │   ├── events.csv
│   │   ├── exon1.bed
│   │   ├── exon2.bed
│   │   ├── exon3.bed
│   │   ├── exon4.bed
│   │   └── metadata.csv
│   └── se
│       ├── events.csv
│       ├── exon1.bed
│       ├── exon2.bed
│       ├── exon3.bed
│       └── metadata.csv
├── junctions
│   ├── metadata.csv
│   └── reads.csv
└── psi
    ├── mxe
    │   └── psi.csv
    ├── outrigger_psi.csv
    └── se
        └── psi.csv

8 directories, 21 files
```

## Check that the found exons are real

Because `outrigger` trusts you, the user, to provide high-quality data, it uses
all available data, including multimapping reads. As a result, many false
positives are expected when detecting novel exons. The best method
for detecting these exons is checking the splice sites using `bedtools` to get
the splice site sequences.

The reason that this step is separate from `outrigger`
is because I want `outrigger` to do one thing (find novel splicing events),
and do that one thing well, rather than doing everything, but doing it poorly.
Plus this requires input of several other files, and IMHO complicates the inputs
to `outrigger`, as it was my goal to make a very simple program.

You'll first want to sort the  `.bed` files. Here's
an example of all the steps you would take for the SE events. The MXE events
are the same, except the folder is `mxe` and you would look at both `exon2` and
`exon3`.

```
cd outrigger_output/index/se/
bedtools sort -i exon2.bed > exon2_sorted.bed
```

Get the upstream flanking sites using `bedtools flank` two nucleotides upstream
(aka to the "left" of the exon) using `-l 2` and `-s` for strand specificity.
You'll want to use `-r 0` to specify no nucleotides to the right.

I obtained the file `~/genomes/mm10/mm10.chrom.sizes`
using the program [`fetchChromSizes`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes)
from [`kentUtils`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/), with
the command `fetchChromSizes mm10 > mm10.chrom.sizes`.

```
bedtools flank -l 2 -r 0 -i exon2_sorted.bed -s -g ~/genomes/mm10/mm10.chrom.sizes > exon2_sorted_2nt_upstream.bed
```

Do the same thing for downstream splice sites, but swap the `-l` and `-r`:

```
bedtools flank -r 2 -l 0 -i exon_sorted2.bed -s -g ~/genomes/mm10/mm10.chrom.sizes > exon2_sorted_2nt_downstream.bed
```

Now you can get the sequence of the splice sites. You'll need a fasta file of the genome sequence.

```
bedtools getfasta -s -tab -fi ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa -bed exon2_sorted_2nt_upstream.bed -fo exon2_sorted_2nt_upstream_sequences.txt
bedtools getfasta -s -tab -fi ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa -bed exon2_sorted_2nt_downstream.bed -fo exon2_sorted_2nt_downstream_sequences.txt
```

You'll get a file like this:

```
$ head exon2_sorted*txt
==> exon2_sorted_2nt_downstream_sequences.txt <==
chr10:128491287-128491289(-)    GT
chr10:128491717-128491719(-)    GT
chr10:128493331-128493333(-)    GT
chr10:128493331-128493333(-)    GT
chr2:136734661-136734663(+) GT
chr2:136734661-136734663(+) GT
chr2:136734661-136734663(+) GT
chr2:136772690-136772692(+) GT
chr2:136772690-136772692(+) GT
chr2:136772690-136772692(+) GT

==> exon2_sorted_2nt_upstream_sequences.txt <==
chr10:128491347-128491349(-)    AG
chr10:128491764-128491766(-)    AG
chr10:128493353-128493355(-)    AG
chr10:128493353-128493355(-)    AG
chr2:136734581-136734583(+) AG
chr2:136734581-136734583(+) AG
chr2:136734581-136734583(+) AG
chr2:136772654-136772656(+) AG
chr2:136772654-136772656(+) AG
chr2:136772654-136772656(+) AG
```

### Filter on only "real" exons

Since the sequences are in the same order as the bed file, you can glue on
these sequences as columns to the original bed file using `cut` and `paste`
and filter for only GT/AG events (or AT/AC for the minor spliceosome).

Here's a one-liner that will get the

```
paste exon2_sorted.bed \
    <(cut -f 2 exon2_sorted_2nt_upstream_sequences.txt) \
    <(cut -f 2 exon2_sorted_2nt_downstream_sequences.txt) \
    > exon2_sorted_with_splice_sites.txt
```

Then the file `exon2_sorted_with_splice_sites.txt` looks like this:

```
chr10	128491289	128491347	isoform1=junction:chr10:128491034-128492058:-|isoform2=junction:chr10:128491348-128492058:-@novel_exon:chr10:128491290-128491347:-@junction:chr10:128491034-128491289:-	.	-	GT	AG
chr10	128491719	128491764	isoform1=junction:chr10:128491034-128492058:-|isoform2=junction:chr10:128491765-128492058:-@novel_exon:chr10:128491720-128491764:-@junction:chr10:128491034-128491719:-	.	-	GT	AG
chr10	128493333	128493353	isoform1=junction:chr10:128492746-128493538:-|isoform2=junction:chr10:128493354-128493538:-@novel_exon:chr10:128493334-128493353:-@junction:chr10:128492746-128493333:-	.	-	GT	AG
chr10	128493333	128493353	isoform1=junction:chr10:128492746-128493538:-|isoform2=junction:chr10:128493354-128493538:-@novel_exon:chr10:128493334-128493353:-@junction:chr10:128492746-128493333:-	.	-	GT	AG
chr2	136734583	136734661	isoform1=junction:chr2:136713601-136756067:+|isoform2=junction:chr2:136713601-136734583:+@novel_exon:chr2:136734584-136734661:+@junction:chr2:136734662-136756067:+	.	+	GT	AG
chr2	136734583	136734661	isoform1=junction:chr2:136713601-136756067:+|isoform2=junction:chr2:136713601-136734583:+@novel_exon:chr2:136734584-136734661:+@junction:chr2:136734662-136756067:+	.	+	GT	AG
chr2	136734583	136734661	isoform1=junction:chr2:136713601-136756067:+|isoform2=junction:chr2:136713601-136734583:+@novel_exon:chr2:136734584-136734661:+@junction:chr2:136734662-136756067:+	.	+	GT	AG
chr2	136772656	136772690	isoform1=junction:chr2:136770175-136773894:+|isoform2=junction:chr2:136770175-136772656:+@novel_exon:chr2:136772657-136772690:+@junction:chr2:136772691-136773894:+	.	+	GT	AG
chr2	136772656	136772690	isoform1=junction:chr2:136770175-136773894:+|isoform2=junction:chr2:136770175-136772656:+@novel_exon:chr2:136772657-136772690:+@junction:chr2:136772691-136773894:+	.	+	GT	AG
chr2	136772656	136772690	isoform1=junction:chr2:136770175-136773894:+|isoform2=junction:chr2:136770175-136772656:+@novel_exon:chr2:136772657-136772690:+@junction:chr2:136772691-136773894:+	.	+	GT	AG
chr2	136773894	136774020	isoform1=junction:chr2:136770175-136777335:+|isoform2=junction:chr2:136770175-136773894:+@exon:chr2:136773895-136774020:+@junction:chr2:136774021-136777335:+	.	+	GT	AG
```

And you can filter for only GT/AG and AT/AC events with:

```
grep -E '(GT\tAG)|(AT\tAC)' exon2_sorted_with_splice_sites.txt
```

This is a relatively boring example because all of the exons would be retained.

## For Developers

How to run the code with the Python debugger. To run the command line functions
such that when they break, you jump into the `pdb` (Python debugger), here is the code:

```
python -m pdb outrigger/commandline.py index \
--sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* \
    --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf
```

Notice that you replace `outrigger` with `python -m pdb outrigger/commandline.py`,
which is relative to this github directory.
