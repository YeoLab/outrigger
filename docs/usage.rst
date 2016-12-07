=====
Usage
=====

To use ``outrigger``, you must perform at least steps 1 and 3:

1. ``outrigger index`` to read all your junctions in and create a database of alternatively spliced exons.
2. ``outrigger validate`` to filter your splicing events for only those with valid introns. (optional)
3. ``outrigger psi`` to calculate percent spliced in on your splicing events.

|OutriggerOverview|

.. |OutriggerOverview| image:: _static/outrigger_overview-300ppi.png
    :target: _static/outrigger_overview-300ppi.png


Each of these commands deserves its own page of explanation, so look at the links below.

Contents:

.. toctree::
    :maxdepth: 2

    subcommands/outrigger_index.rst
    subcommands/outrigger_validate.rst
    subcommands/outrigger_psi.rst
