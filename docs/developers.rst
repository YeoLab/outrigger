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

If you want to run a smallish example with GENCODE GTF files and a mouse
genome, do:

::

    make tasic2016

To run this with different numbers of parallel processing cores, do:

If you want to run a smallish example with GENCODE GTF files and a mouse
genome, specify with ``N_JOBS=X``, where ``X`` is the number of jobs you
want. By default, this uses ``-1`` jobs, which means to use the maximum
number of processors available.

::

    make tasic2016 N_JOBS=8

If you want to run an example with ENSEMBL GTF files, do:

::

    make arabdopsis

By default, Travis-CI checks for coverage and that the Arabdopsis example runs.

::

    script:
    - make coverage
    - make arabdopsis

Checking code style (linting)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to issues with ``bioconda`` builds not allowing for setuptools-installing
packages, the ``flake8`` packages used to enforce PEP8 code style and practices
is not part of the ``requirements.txt`` or ``environment.yml``. So, instead, on
Travis, we create an environment and recommend for developers to do the same.

From the ``outrigger`` root directory, where there is a ``Makefile`` defining
``make lint``, do:

::

    conda create -n lint-env --yes flake8
    source activate lint-env
    make lint
    deactivate lint-env
