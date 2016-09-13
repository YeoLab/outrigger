.. :changelog:

History
=======

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
