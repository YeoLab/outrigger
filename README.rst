Outrigger
=========

|Outrigger logo|

|Build Status|\ |image2|\ |Coverage Status|

Outrigger is a program which uses junction reads from RNA seq data, and
a graph database to create a *de novo* alternative splicing annotation
with a graph database, and quantify percent spliced-in (Psi) of the
events.

-  Free software: BSD license

Features
--------

-  Finds novel splicing events, including novel exons!
   (``outrigger index``)

   -  Currently only works with ``SJ.out.tab`` files from the
      `STAR <https://github.com/alexdobin/STAR>`__ aligner

-  (optional) Validates that exons have correct splice sites, e.g. GT/AG
   and AT/AC for mammalian systems (``outrigger validate``)
-  Calculate "percent spliced-in" (Psi/Ψ) scores for all your samples
   given the validated events (or the original events if you opted not
   to validate)

Installation
------------

To install ``outrigger``, we recommend using the `Anaconda Python
Distribution <http://anaconda.org/>`__ and creating an environment.

You'll want to add the ```bioconda`` <https://bioconda.github.io/>`__
channel to make installing ```bedtools`` <bedtools.readthedocs.io>`__
and its Python wrapper,
```pybedtools`` <https://daler.github.io/pybedtools/>`__ easy.

::

    conda config --add channels r
    conda config --add channels bioconda

Create an environment called ``outrigger-env``. Python 2.7, Python 3.4,
and Python 3.5 are supported.

::

    conda create -n outrigger-env pandas pybedtools gffutils biopython bedtools joblib

Now activate that environment using ``source activate outrigger-env``
and install ``outrigger`` from PyPI, using ``pip``:

::

    source activate outrigger-env
    pip install outrigger

To check that it installed properly, try the command with the help
option (``-h``), ``outrigger -h``. The output should look like this:

::

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

Bleeding edge code from Github (here)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For advanced users, if you have `git <https://git-scm.com/>`__ and
`Anaconda Python <https://www.continuum.io/downloads>`__ installed, you
can:

#. Clone this repository
#. Change into that directory
#. Create an environment with the necessary packages from Anaconda
#. Activate the environment
#. Install remaining packages from PyPI
   (```graphlite`` <https://github.com/eugene-eeo/graphlite>`__ is only
   available on PyPI, not as a ``conda`` package)
#. Install this package

These steps are shown in code below.

::

    git clone git@github.com:YeoLab/outrigger
    cd outrigger
    conda create --name outrigger --yes --file conda_requirements.txt --channel bioconda
    source activate outrigger
    pip install -r requirements.txt
    pip install .

Quick start
-----------

If you just want to know how to run this on your data with the default
parameters, start here. Let's say you performed your alignment in the
folder called ``~/projects/tasic2016/analysis/tasic2016_v1``, and that's
where your ``SJ.out.tab`` files from the STAR aligner are (they're
output into the same folder as the ``.bam`` files). First you'll need to
change directories to that folder with ``cd``.

::

    cd ~/projects/tasic2016/analysis/tasic2016_v1

Then you need find all alternative splicing events, which you do by
running ``outrigger index`` on the splice junction files and the gtf.
Here is an example command:

::

    outrigger index --sj-out-tab *SJ.out.tab \
        --gtf /projects/ps-yeolab/genomes/mm10/gencode/m10/gencode.vM10.annotation.gtf

Next, you'll want to validate that the splicing events you found follow
biological rules, such as being containing GT/AG (mammalian major
spliceosome) or AT/AC (mammalian minor splicesome) sequences. To do
that, you'll need to provide the genome name (e.g. ``mm10``) and the
genome sequences. An example command is below:

::

    outrigger validate --genome mm10 \
        --fasta /projects/ps-yeolab/genomes/mm10/GRCm38.primary_assembly.genome.fa

Finally, you can calculate percent spliced in (Psi) of your splicing
events! Thankfully this is very easy:

::

    outrigger psi

It should be noted that ALL of these commands should be performed in the
same directory, so no moving.

Quick start summary
~~~~~~~~~~~~~~~~~~~

Here is a summary the commands in the order you would use them for
outrigger!

::

    cd ~/projects/tasic2016/analysis/tasic2016_v1
    outrigger index --sj-out-tab *SJ.out.tab \
        --gtf /projects/ps-yeolab/genomes/mm10/gencode/m10/gencode.vM10.annotation.gtf
    outrigger validate --genome mm10 \
        --fasta /projects/ps-yeolab/genomes/mm10/GRCm38.primary_assembly.genome.fa
    outrigger psi

This will create a folder called ``outrigger_output``, which at the end
should look like this:

::

    $ tree outrigger_output
    outrigger_output
    ├── index
    │   ├── gtf
    │   │   ├── gencode.vM10.annotation.gtf
    │   │   ├── gencode.vM10.annotation.gtf.db
    │   │   └── novel_exons.gtf
    │   ├── junction_exon_direction_triples.csv
    │   ├── mxe
    │   │   ├── events.csv
    │   │   ├── exon1.bed
    │   │   ├── exon2.bed
    │   │   ├── exon3.bed
    │   │   ├── exon4.bed
    │   │   ├── splice_sites.csv
    │   │   └── validated
    │   │       └── events.csv
    │   └── se
    │       ├── events.csv
    │       ├── exon1.bed
    │       ├── exon2.bed
    │       ├── exon3.bed
    │       ├── splice_sites.csv
    │       └── validated
    │           └── events.csv
    ├── junctions
    │   ├── metadata.csv
    │   └── reads.csv
    └── psi
        ├── mxe
        │   └── psi.csv
        ├── outrigger_psi.csv
        └── se
            └── psi.csv

    10 directories, 22 files

Commands
--------

Here's an in-depth look at the commands of \`outrigger.

``index``: Build a *de novo* splicing annotation index of events custom to *your* data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The "help" output of the two programs tries to be explicit about what is
required to run ``outrigger``. Below is the output of when you use the
command, ``outrigger index --help``

::

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

Example ``outrigger index`` command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Included in this repository is a subset of the 1809 cells from `"Adult
mouse cortical cell taxonomy revealed by single cell transcriptomics."
by Tasic et al, Nature Neuroscience
(2016) <http://www.ncbi.nlm.nih.gov/pubmed/26727548>`__. There splice
junction output files from the `STAR
aligner <https://github.com/alexdobin/STAR>`__ from the 43
"``CAV_LP_Ipsi_tdTpos``\ " cells, plus a subset of the `GENCODE
M10 <http://www.gencodegenes.org/mouse_releases/10.html>`__ (Version M10
(January 2016 freeze, GRCm38) - Ensembl 85) mouse annotation.

To run this program with the included example data, from the
``outrigger`` directory where you cloned ``outrigger`` (this is
important because the locations of the files is relative to that
directory), run this command:

::

    outrigger index \
        --sj-out-tab outrigger/tests/data/tasic2016/unprocessed/sj_out_tab/* \
        --gtf outrigger/tests/data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf

*Note: the backslashes (``\``, like a tree that's falling backwards
relative to right-to-left reading) to tell the computer that you're not
done telling it what to do, so the line continues, and to help the code
be a little more human-readable. The above code is read by the computer
exactly the same as the one-liner below, but is easier for us humans to
read.*

::

    outrigger index --sj-out-tab outrigger/tests/data/tasic2016/unprocessed/sj_out_tab/* --gtf outrigger/tests/data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf

This is equivalent to the below command, which specifies all the other
arguments with the default values.

::

    outrigger index \
        --sj-out-tab outrigger/tests/data/tasic2016/unprocessed/sj_out_tab/* \
        --gtf outrigger/tests/data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf \
        --output ./outrigger_output --min-reads 10

The output of this command is:

::

    $ outrigger index --sj-out-tab example_data/tasic2016/unprocessed/sj_out_tab/* --gtf example_data/tasic2016/unprocessed/gtf/snap25_myl6.gtf
    2016-08-12 11:24:03 Reading SJ.out.files and creating a big splice junction table of reads spanning exon-exon junctions...
    2016-08-12 11:24:03 Writing ./outrigger_output/junction_reads/reads.csv ...

    2016-08-12 11:24:03     Done.
    2016-08-12 11:24:03 Creating splice junction metadata of merely where junctions start and stop
    2016-08-12 11:24:03     Done.
    2016-08-12 11:24:03 Getting junction-direction-exon triples for graph database ...
    2016-08-12 11:24:03 Starting annotation of all junctions with known neighboring exons ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Writing junction-exon-direction triples to ./outrigger_output/index/junction_exon_direction_triples.csv...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Populating graph database of the junction-direction-exon triples ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Finding all skipped exon (SE) events ...
    2016-08-12 11:24:04 Trying out 25 exons ...
    2016-08-12 11:24:04     1/25 exons tested (4.0%)
    2016-08-12 11:24:04     2/25 exons tested (8.0%)
    2016-08-12 11:24:04     3/25 exons tested (12.0%)
    2016-08-12 11:24:04     4/25 exons tested (16.0%)
    2016-08-12 11:24:04     5/25 exons tested (20.0%)
    2016-08-12 11:24:04     6/25 exons tested (24.0%)
    2016-08-12 11:24:04     7/25 exons tested (28.0%)
    2016-08-12 11:24:04     8/25 exons tested (32.0%)
    2016-08-12 11:24:04     9/25 exons tested (36.0%)
    2016-08-12 11:24:04     10/25 exons tested (40.0%)
    2016-08-12 11:24:04     11/25 exons tested (44.0%)
    2016-08-12 11:24:04     12/25 exons tested (48.0%)
    2016-08-12 11:24:04     13/25 exons tested (52.0%)
    2016-08-12 11:24:04     14/25 exons tested (56.0%)
    2016-08-12 11:24:04     15/25 exons tested (60.0%)
    2016-08-12 11:24:04     16/25 exons tested (64.0%)
    2016-08-12 11:24:04     17/25 exons tested (68.0%)
    2016-08-12 11:24:04     18/25 exons tested (72.0%)
    2016-08-12 11:24:04     19/25 exons tested (76.0%)
    2016-08-12 11:24:04     20/25 exons tested (80.0%)
    2016-08-12 11:24:04     21/25 exons tested (84.0%)
    2016-08-12 11:24:04     22/25 exons tested (88.0%)
    2016-08-12 11:24:04     23/25 exons tested (92.0%)
    2016-08-12 11:24:04     24/25 exons tested (96.0%)
    2016-08-12 11:24:04     25/25 exons tested (100.0%)
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Writing 1 SE events to ./outrigger_output/index/se/events.csv ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Making metadata file of SE events, annotating them with GTF attributes ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Getting exon and intron lengths of alternative events ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Combining lengths and attributes into one big dataframe ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Writing SE metadata to ./outrigger_output/index/se/metadata.csv ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Finding all mutually exclusive exon (MXE) events ...
    2016-08-12 11:24:04 Trying out 25 exons ...
    2016-08-12 11:24:04     1/25 exons tested (4.0%)
    2016-08-12 11:24:04     2/25 exons tested (8.0%)
    2016-08-12 11:24:04     3/25 exons tested (12.0%)
    2016-08-12 11:24:04     4/25 exons tested (16.0%)
    2016-08-12 11:24:04     5/25 exons tested (20.0%)
    2016-08-12 11:24:04     6/25 exons tested (24.0%)
    2016-08-12 11:24:04     7/25 exons tested (28.0%)
    2016-08-12 11:24:04     8/25 exons tested (32.0%)
    2016-08-12 11:24:04     9/25 exons tested (36.0%)
    2016-08-12 11:24:04     10/25 exons tested (40.0%)
    2016-08-12 11:24:04     11/25 exons tested (44.0%)
    2016-08-12 11:24:04     12/25 exons tested (48.0%)
    2016-08-12 11:24:04     13/25 exons tested (52.0%)
    2016-08-12 11:24:04     14/25 exons tested (56.0%)
    2016-08-12 11:24:04     15/25 exons tested (60.0%)
    2016-08-12 11:24:04     16/25 exons tested (64.0%)
    2016-08-12 11:24:04     17/25 exons tested (68.0%)
    2016-08-12 11:24:04     18/25 exons tested (72.0%)
    2016-08-12 11:24:04     19/25 exons tested (76.0%)
    2016-08-12 11:24:04     20/25 exons tested (80.0%)
    2016-08-12 11:24:04     21/25 exons tested (84.0%)
    2016-08-12 11:24:04     22/25 exons tested (88.0%)
    2016-08-12 11:24:04     23/25 exons tested (92.0%)
    2016-08-12 11:24:04     24/25 exons tested (96.0%)
    2016-08-12 11:24:04     25/25 exons tested (100.0%)
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Writing 1 MXE events to ./outrigger_output/index/mxe/events.csv ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Making metadata file of MXE events, annotating them with GTF attributes ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Getting exon and intron lengths of alternative events ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Combining lengths and attributes into one big dataframe ...
    2016-08-12 11:24:04     Done.
    2016-08-12 11:24:04 Writing MXE metadata to ./outrigger_output/index/mxe/metadata.csv ...
    2016-08-12 11:24:04     Done.

``outrigger_index`` Outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The above commands will create a folder called ``outrigger_index`` in
the folder you ran the command from, with the following structure

::

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

``validate``: Check that the found exons are real
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example command assumes that you have a ``mm10`` genome fasta file
located at
``~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa`` and a
chromosome sizes file located at ``~/genomes/mm10/mm10.chrom.sizes``

::

    outrigger validate -f ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa -g ~/genomes/mm10/mm10.chrom.sizes

``psi``: Calculate percent spliced-in (Psi/Ψ) scores for your data from the splicing events you created
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

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

::

    outrigger psi

The above command is equivalent to specifying all the arguments with
their default values:

::

    outrigger psi --index ./outrigger_index --min-reads 10

``outrigger_psi`` Outputs
^^^^^^^^^^^^^^^^^^^^^^^^^

Now the ``outrigger_output`` folder has ``psi`` subfolder, with the MXE
and SE events separate.

::

    $ tree outrigger_output
    outrigger_output
    ├── index
    │   ├── gtf
    │   │   ├── gencode.vM10.annotation.subset.gtf
    │   │   ├── gencode.vM10.annotation.subset.gtf.db
    │   │   └── novel_exons.gtf
    │   ├── junction_exon_direction_triples.csv
    │   ├── mxe
    │   │   ├── events.csv
    │   │   ├── exon1.bed
    │   │   ├── exon2.bed
    │   │   ├── exon3.bed
    │   │   ├── exon4.bed
    │   │   ├── splice_sites.csv
    │   │   └── validated
    │   │       └── events.csv
    │   └── se
    │       ├── events.csv
    │       ├── exon1.bed
    │       ├── exon2.bed
    │       ├── exon3.bed
    │       ├── splice_sites.csv
    │       └── validated
    │           └── events.csv
    ├── junctions
    │   ├── metadata.csv
    │   └── reads.csv
    └── psi
        ├── mxe
        │   └── psi.csv
        ├── outrigger_psi.csv
        └── se
            └── psi.csv

    10 directories, 22 files

For Developers
--------------

How to run the code with the Python debugger. To run the command line
functions such that when they break, you jump into the ``pdb`` (Python
debugger), here is the code:

::

    python -m pdb outrigger/commandline.py index \
    --sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* \
        --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf

Notice that you replace ``outrigger`` with
``python -m pdb outrigger/commandline.py``, which is relative to this
github directory.

.. |Outrigger logo| image:: https://raw.githubusercontent.com/YeoLab/outrigger/master/logo/logo_v1.png
.. |Build Status| image:: https://travis-ci.org/YeoLab/outrigger.svg?branch=master
   :target: https://travis-ci.org/YeoLab/outrigger
.. |image2| image:: https://img.shields.io/pypi/v/outrigger.svg
   :target: https://pypi.python.org/pypi/outrigger
.. |Coverage Status| image:: https://coveralls.io/repos/YeoLab/outrigger/badge.svg?branch=master&service=github
   :target: https://coveralls.io/github/YeoLab/outrigger?branch=master
