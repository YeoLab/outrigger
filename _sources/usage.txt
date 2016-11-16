=====
Usage
=====

To use ``outrigger``, you must perform at least steps 1 and 3:

1. ``outrigger index`` to read all your junctions in and create a splicing database.
2. ``outrigger validate`` to filter your splicing events for only those with valid introns. (optional)
3. ``outrigger psi`` to calculate percent spliced in on your splicing events.

.. include:: subcommands/outrigger_index
.. include:: subcommands/outrigger_validate
.. include:: subcommands/outrigger_psi
