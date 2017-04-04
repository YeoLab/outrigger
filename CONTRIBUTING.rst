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

Tips for writing documentation
------------------------------

To keep documentation writing separate from code and test writing, I have an
environment only for documentation writing that was created in this way:

```
# Go to the folder where you cloned outrigger
cd outrigger
conda create -n outrigger-documentation docs/requirements.txt
```

Now activate the environment and build the docs. This command will both build
what you have and open a web browser so you can see the latest version.

```
source activate outrigger-documentation
make docs
```

Make changes as you want in the `docs` folder. Once you're ready to deploy
them to the `yeolab.github.io/outrigger` website, you'll need to change to the
`docs` directory and use the `make deploy` directive there to build the docs
and push the HTML files to the `gh-pages` branch. Here is an example with
output:

```
(outrigger-documentation) ➜  docs git:(documentation) ✗ cd docs
(outrigger-documentation) ➜  docs git:(documentation) ✗ make deploy
echo "Preparing gh_pages deployment..."
Preparing gh_pages deployment...
echo "Pulling any updates from Github Pages..."
Pulling any updates from Github Pages...
cd _deploy; git pull;
remote: Counting objects: 2039, done.
remote: Compressing objects: 100% (127/127), done.
remote: Total 2039 (delta 462), reused 402 (delta 401), pack-reused 1510
Receiving objects: 100% (2039/2039), 164.78 MiB | 9.19 MiB/s, done.
Resolving deltas: 100% (1206/1206), completed with 84 local objects.
From github.com:YeoLab/outrigger
 * [new branch]      checklist  -> origin/checklist
 + c68a887...04316b6 documentation -> origin/documentation  (forced update)
 * [new branch]      flake8-weirdness -> origin/flake8-weirdness
 * [new branch]      logo       -> origin/logo
 * [new branch]      low_memory -> origin/low_memory
   82b6296..9c80211  master     -> origin/master
 * [new branch]      merge_exon_finding_and_junction_adjacencies -> origin/merge_exon_finding_and_junction_adjacencies
 * [new branch]      py3.6      -> origin/py3.6
 * [new branch]      segfault   -> origin/segfault
 * [new branch]      strandedness -> origin/strandedness
 * [new branch]      v1.0.0rc1  -> origin/v1.0.0rc1
 * [new tag]         v1.0.0     -> v1.0.0
Already up-to-date.
mkdir -p _deploy/
echo "Copying files from '_build/html/.' to '_deploy/'"
Copying files from '_build/html/.' to '_deploy/'
cp -r _build/html/. _deploy/
echo "Deploying on github pages now..."
Deploying on github pages now...
cd _deploy; git add -A; git commit -m "docs updated at `date -u`";\
		git push origin  --quiet
[gh-pages fa0420c] docs updated at Tue Apr  4 01:31:17 UTC 2017
 134 files changed, 103905 insertions(+), 2704 deletions(-)
 create mode 100644 _images/outrigger_index-300ppi.png
 create mode 100644 _images/outrigger_overview-300ppi.png
 create mode 100644 _images/outrigger_psi-1x.png
 create mode 100644 _images/outrigger_validate-1x.png
 rewrite _modules/outrigger/psi/compute.html (68%)
 create mode 100644 _sources/license.txt
 rewrite _sources/releases/v0.2.7.txt (83%)
 create mode 100644 _sources/releases/v1.0.0.txt
 create mode 100644 _sources/releases/v1.0.1.txt
 create mode 100644 _static/exon_vs_junction_start_stop.ai
 create mode 100644 _static/logo-150ppi.png
 create mode 100644 _static/logo-1x.png
 create mode 100644 _static/logo-300ppi.png
 create mode 100644 _static/logo.ai
 create mode 100644 _static/logo.svg
 create mode 100644 _static/outrigger_franken-events.ai
 create mode 100644 _static/outrigger_index-150ppi.png
 create mode 100644 _static/outrigger_index-1x.png
 create mode 100644 _static/outrigger_index-300ppi.png
 create mode 100644 _static/outrigger_index.ai
 create mode 100644 _static/outrigger_index.svg
 create mode 100644 _static/outrigger_index_flanking_exons.ai
 create mode 100644 _static/outrigger_index_mxe.ai
 create mode 100644 _static/outrigger_index_se.ai
 create mode 100644 _static/outrigger_overview-150ppi.png
 create mode 100644 _static/outrigger_overview-1x.png
 create mode 100644 _static/outrigger_overview-300ppi.png
 create mode 100644 _static/outrigger_overview.svg
 create mode 100644 _static/outrigger_overview_v1.ai
 create mode 100644 _static/outrigger_overview_v1.svg
 create mode 100644 _static/outrigger_overview_v2.ai
 create mode 100644 _static/outrigger_psi-150ppi.png
 rewrite _static/outrigger_psi-1x.png (98%)
 rewrite _static/outrigger_psi-300ppi.png (92%)
 create mode 100644 _static/outrigger_psi_pathological_cases.ai
 create mode 100644 _static/outrigger_psi_v1.ai
 create mode 100644 _static/outrigger_psi_v2.ai
 create mode 100644 _static/outrigger_psi_v3.ai
 create mode 100644 _static/outrigger_validate-150ppi.png
 create mode 100644 _static/outrigger_validate-1x.png
 create mode 100644 _static/outrigger_validate-300ppi.png
 create mode 100644 _static/outrigger_validate.ai
 create mode 100644 _static/outrigger_validate.svg
 create mode 100644 license.html
 rewrite objects.inv (97%)
 create mode 100644 releases/v1.0.0.html
 create mode 100644 releases/v1.0.1.html
 rewrite searchindex.js (99%)
echo "Github Pages deploy was completed at `date -u`"
Github Pages deploy was completed at Tue Apr  4 01:31:24 UTC 2017```

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
