.. :changelog:

History
=======

v0.2.8 (October 23rd, 2016)
---------------------------

Updated README/HISTORY files


v0.2.7 (October 23rd, 2016)
---------------------------

New features
~~~~~~~~~~~~

- Added ``outrigger validate`` command to check for canonical splice sites
  by default: ``GT/AG`` (U1, major spliceosome) and ``AT/AC``
  (U12, minor spliceosome). Both of these are user-adjustable as they are only
  the standard for mammalian genomes.

API changes
~~~~~~~~~~~

- Added ``--resume`` and ``--force`` options to ``outrigger index`` to prevent
  the overwriting of interrupted indexing operations, or to force overwriting.
  By default, ``outrigger`` complains and cowardly exits.

Bug fixes
~~~~~~~~~

- Support ENSEMBL gtf files which specify chromsome names with a number, e.g.
  ``4`` instead of ``chr4``. Thank you to lcscs12345_ for pointing this out!

Miscellaneous
~~~~~~~~~~~~~

- Added version info with ``outrigger --version``
- Sped up gffutils queries and event finding by running ``ANALYZE`` on SQLite
  databases.


.. _lcscs12345: https://github.com/lcscs12345


v0.2.6 (September 15th, 2016)
-----------------------------

This is a non-breaking patch release

Bug fixes
~~~~~~~~~

- Wasn't concatenating exons properly after parallelizing


Miscellaneous
~~~~~~~~~~~~~

- Clarified ``.gtf`` file example for directory output



v0.2.5 (September 14th, 2016)
-----------------------------


Bug fixes
~~~~~~~~~

- Added ``joblib`` to requirements


v0.2.4 (September 14th, 2016)
-----------------------------

This is a non-breaking patch release of ``outrigger``.

New features
~~~~~~~~~~~~

- **Actually** parallelized exon finding for novel exons. Before had written the code and tested the non-parallelized version but now using actually parallelized version!


Bug fixes
~~~~~~~~~

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
