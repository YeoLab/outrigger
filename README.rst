Outrigger
=========

|OutriggerLogo|

|BuildStatus|\ |codecov|\ |PyPI
versions|\ |PythonVersionCompatibility|\ |license|

Outrigger is a program which uses junction reads from RNA seq data, and
a graph database to create a *de novo* alternative splicing annotation
with a graph database, and quantify percent spliced-in (Psi) of the
events.

-  Free software: BSD license

Features
--------

-  Finds novel splicing events, including novel exons!
   (``outrigger index``) from ``.bam`` files
-  (optional) Validates that exons have correct splice sites, e.g. GT/AG
   and AT/AC for mammalian systems (``outrigger validate``)
-  Calculate "percent spliced-in" (Psi/Ψ) scores for all your samples
   given the validated events (or the original events if you opted not
   to validate) via ``outrigger psi``

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

For Developers
--------------

How to run with the Python debugger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

How to run the tests
~~~~~~~~~~~~~~~~~~~~

| If you want to run the tests without calculating what percentage of
| lines are covered in the test suite, run:

::

    make test

| If you want to run the tests and see which lines are covered by tests
and get
| an overall percentage of test coverage, run:

::

    make coverage

If you want to run an example with ENSEMBL GTF files, do:

::

    make arabdopsis

By default, Travis-CI does all three:

::

    script:
    - make coverage
    - make lint
    - make arabdopsis

.. |OutriggerLogo| image:: https://raw.githubusercontent.com/YeoLab/outrigger/master/logo/logo_v1.png
.. |BuildStatus| image:: https://travis-ci.org/YeoLab/outrigger.svg?branch=master
   :target: https://travis-ci.org/YeoLab/outrigger
.. |codecov| image:: https://codecov.io/gh/YeoLab/outrigger/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/YeoLab/outrigger
.. |PyPI versions| image:: https://img.shields.io/pypi/v/outrigger.svg
   :target: https://pypi.python.org/pypi/outrigger
.. |PythonVersionCompatibility| image:: https://img.shields.io/pypi/pyversions/outrigger.svg
   :target: https://pypi.python.org/pypi/outrigger
.. |license| image:: https://img.shields.io/github/pypi/l/outrigger.svg
   :target: 
