.. :changelog:

History
=======

v1.1.0 (June 28th, 2017)
------------------------

This is a minor release to ```outrigger``.

Bug fixes
~~~~~~~~~

- Fixed `UNIQUE ID` error that happened somewhat stochastically when adding new exons to the database


Miscellaneous
~~~~~~~~~~~~~

- Explicitly added Python 3.6 compatibility
- Change logo location to `logo/` instead of `logo/v1` since there's only one
  version anyway...


v1.0.0 (April 3rd, 2017)
------------------------

This is the first major release of ``outrigger``!!!

v1.0.0 New features
~~~~~~~~~~~~~~~~~~~

- Parallelized event across chromosomes
- Added ``--low-memory`` flag for ``index``, ``validate``, and ``psi`` commands
  to use a smaller memory footprint when reading CSV files.
- Added ``--splice-types`` option to specify only one kind of splicing you'd
  like to find
- So the user can double-check the Psi calculation, create a ``summary.csv``
  file indicating the number of reads found at each junction, for all samples
  - This also shows which "Case" corresponds to each event in each sample, so you can see whether there were sufficient or insufficient reads on the junctions of each event, and how ``outrigger`` judged it.
- Added functions to extract constitutive and alternative exons separately

v1.0.0 Bug fixes
~~~~~~~~~~~~~~~~

- Fixed a bug that stalled on ``.bam`` files while counting the junctions

v1.0.0 Miscellaneous
~~~~~~~~~~~~~~~~~~~~

- Added ``GC/AG`` to valid splice sites


v0.2.9 (November 11th, 2016)
----------------------------

This is a non-breaking release with many speed improvements, and upgrade is
recommended.


v0.2.9 New features
~~~~~~~~~~~~~~~~~~~

- Add ``bam`` alignment files as input option


Miscellaneous
~~~~~~~~~~~~~

- Parallelized Psi calculation, the exact number of processors can be specified
  with ``--n-jobs``, and by default, ``--n-jobs`` is ``-1``, which means use as
  many processors as are available.


v0.2.8 (October 23rd, 2016)
---------------------------

Updated README/HISTORY files


v0.2.7 (October 23rd, 2016)
---------------------------

v0.2.7 New features
~~~~~~~~~~~~~~~~~~~

- Added ``outrigger validate`` command to check for canonical splice sites
  by default: ``GT/AG`` (U1, major spliceosome) and ``AT/AC``
  (U12, minor spliceosome). Both of these are user-adjustable as they are only
  the standard for mammalian genomes.

v0.2.7 API changes
~~~~~~~~~~~~~~~~~~

- Added ``--resume`` and ``--force`` options to ``outrigger index`` to prevent
  the overwriting of interrupted indexing operations, or to force overwriting.
  By default, ``outrigger`` complains and cowardly exits.

v0.2.7 Bug fixes
~~~~~~~~~~~~~~~~

- Support ENSEMBL gtf files which specify chromsome names with a number, e.g.
  ``4`` instead of ``chr4``. Thank you to lcscs12345_ for pointing this out!

v0.2.7 Miscellaneous
~~~~~~~~~~~~~~~~~~~~

- Added version info with ``outrigger --version``
- Sped up gffutils queries and event finding by running ``ANALYZE`` on SQLite
  databases.


.. _lcscs12345: https://github.com/lcscs12345


v0.2.6 (September 15th, 2016)
-----------------------------

This is a non-breaking patch release

v0.2.6 Bug fixes
~~~~~~~~~~~~~~~~

- Wasn't concatenating exons properly after parallelizing


v0.2.6 Miscellaneous
~~~~~~~~~~~~~~~~~~~~

- Clarified ``.gtf`` file example for directory output



v0.2.5 (September 14th, 2016)
-----------------------------


v0.2.5 Bug fixes
~~~~~~~~~~~~~~~~

- Added ``joblib`` to requirements


v0.2.4 (September 14th, 2016)
-----------------------------

This is a non-breaking patch release of ``outrigger``.

v0.2.4 New features
~~~~~~~~~~~~~~~~~~~

- **Actually** parallelized exon finding for novel exons. Before had written the code and tested the non-parallelized version but now using actually parallelized version!


v0.2.4 Bug fixes
~~~~~~~~~~~~~~~~

- Don't need to turn on ``--debug`` command for outrigger to even run



v0.2.3 (September 13th, 2016)
-----------------------------

This is a patch release of outrigger, with non-breaking changes from the
previous one.


Bug fixes
~~~~~~~~~

- Subfolders get copied when installing
- Add test for checking that ``outrigger -h`` command works


v0.2.2 (September 12th, 2016)
-----------------------------

This is a point release which includes the ``index`` submodule in the ``__all__`` statement.


v0.2.1 (September 12th, 2016)
-----------------------------

This is a point release which actually includes the ``requirements.txt`` file that specifies which packages ``outrigger`` depends on.


v0.2.0 (September 9th, 2016)
----------------------------

This is the second release of ``outrigger``!

New features
~~~~~~~~~~~~

- Parallelized exon finding for novel exons
- Added ``outrigger validate`` command to check that your new exons have proper splice sites (e.g. GT/AG and AT/AC)
- Added more test data for other event types (even though we don't detect them yet)


v0.1.0 (May 25, 2016)
---------------------

This is the initial release of ``outrigger``
