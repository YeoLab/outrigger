============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/YeoLab/outrigger/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Outrigger could always use more documentation, whether as part of the
official Outrigger docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/YeoLab/outrigger/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up ``outrigger`` for local development.

1. Fork the ``outrigger`` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/outrigger.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv outrigger
    $ cd outrigger/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox::

    $ flake8 outrigger tests
    $ py.test
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.


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



Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.6, 2.7, 3.3, and 3.4, and for PyPy. Check
   https://travis-ci.org/olgabot/outrigger/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

    $ python -m unittest tests.test_outrigger
