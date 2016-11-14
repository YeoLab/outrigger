``validate``: Check that the found exons are real
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example command assumes that you have a ``mm10`` genome fasta file
located at
``~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa`` and a
chromosome sizes file located at ``~/genomes/mm10/mm10.chrom.sizes``

::

    outrigger validate -f ~/genomes/mm10/gencode/m10/GRCm38.primary_assembly.genome.fa -g ~/genomes/mm10/mm10.chrom.sizes

