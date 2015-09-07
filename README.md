# PoshSplice

[![](https://img.shields.io/travis/olgabot/poshsplice.svg)](https://travis-ci.org/olgabot/poshsplice)[![](https://img.shields.io/pypi/v/poshsplice.svg)](https://pypi.python.org/pypi/poshsplice)

Exon ontology is an attempt to catalogue scripts used to annotate exons and their functions

* Free software: BSD license
* Documentation: https://olgabot.github.io/poshsplice

## Features

### Aggregate splice junctions into splicing events

Using just junction information, we can annotate these kinds of splicing 
events:

- Alternative first exon (AFE)
- Skipped exon (SE)
- Alternative 5' splice site (A5SS)
- Alternative 3' splice site (A3SS)
- Mutually exclusive exon (MXE)
- Twin cassette
- Alternative last exon

For a demonstration, see the transcripts below.

![](docs/test_transcripts.png)

### Read HMMScan output

After running HMMScan, the output is extensive and helpful, but non-trivial to parse into a dataframe-style format
since the separators are variable lengths of spaces, and the domain descriptions also have spaces, so parsing is
ambiguous. Here's an example output file, of using the [RBFOX2](http://en.wikipedia.org/wiki/RBM9) protein sequence
as a query, and scanning against the Pfam-A (manually curated) families from the Pfam version 27:

```
#                                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name            accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- -----  -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
Fox-1_C              PF12414.3     93 sp|O43251|RFOX2_HUMAN -            390   3.2e-39  133.2  29.5   1   2      0.23   6.7e+02    0.7   0.0    14    48   177   213   166   243 0.66 Calcitonin gene-related peptide regulator C terminal
Fox-1_C              PF12414.3     93 sp|O43251|RFOX2_HUMAN -            390   3.2e-39  133.2  29.5   2   2   8.9e-42   2.6e-38  130.2  27.3     2    93   265   362   264   362 0.97 Calcitonin gene-related peptide regulator C terminal
RRM_1                PF00076.17    70 sp|O43251|RFOX2_HUMAN -            390     8e-19   67.0   0.1   1   1   5.9e-22   1.7e-18   65.9   0.1     2    70   124   191   123   191 0.97 RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)
RRM_6                PF14259.1     70 sp|O43251|RFOX2_HUMAN -            390   2.4e-15   56.2   0.1   1   1   1.4e-18   4.3e-15   55.4   0.1     1    70   123   191   123   191 0.95 RNA recognition motif (a.k.a. RRM, RBD, or RNP domain)
RRM_5                PF13893.1     56 sp|O43251|RFOX2_HUMAN -            390   8.1e-11   41.6   0.1   1   1   5.9e-14   1.8e-10   40.5   0.1     1    54   137   193   137   195 0.90 RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)
RRM_3                PF08777.6    105 sp|O43251|RFOX2_HUMAN -            390     0.084   12.7   0.0   1   1   6.7e-05       0.2   11.5   0.0    17    79   136   202   127   206 0.83 RNA binding motif
#
# Program:         hmmscan
# Version:         3.1b1 (May 2013)
# Pipeline mode:   SCAN
# Query file:      /projects/ps-yeolab/genomes/pfam/release_27/ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam27.0/RBFOX2_human.fasta
# Target file:     /projects/ps-yeolab/genomes/pfam/release_27/ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm
# Option settings: hmmscan --domtblout RBFOX2_human_pfam.txt --noali --notextw /projects/ps-yeolab/genomes/pfam/release_27/ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm /projects/ps-yeolab/genomes/pfam/release_27/ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam27.0/RBFOX2_human.fasta
# Current dir:     /home/obotvinnik/ipython_notebook/miso
# Date:            Tue Jan 27 18:56:23 2015
# [ok]
```

You can read in this file using `read_hmmscan_out`, as such:

```python
from poshsplice.hmmscan import read_hmmscan
hmmscan_df = read_hmmscan('rbfox2_hmmscan.txt')
```

This will give you a pandas dataframe which looks like this:

<table border="1" class="dataframe">\n  <thead>\n    <tr style="text-align: right;">\n      <th></th>\n      <th>target_name</th>\n      <th>target_accession</th>\n      <th>target_length</th>\n      <th>query_name</th>\n      <th>query_accession</th>\n      <th>query_length</th>\n      <th>sequence_e_value</th>\n      <th>sequence_score</th>\n      <th>sequence_bias</th>\n      <th>domain_number</th>\n      <th>domain_total</th>\n      <th>domain_conditional_e_value</th>\n      <th>domain_independent_e_value</th>\n      <th>domain_score</th>\n      <th>domain_bias</th>\n      <th>target_start</th>\n      <th>target_stop</th>\n      <th>query_start</th>\n      <th>query_stop</th>\n      <th>query_domain_envelope_start</th>\n      <th>query_domain_envelope_stop</th>\n      <th>mean_posterior_probability</th>\n      <th>target_description</th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr>\n      <th>0</th>\n      <td>Fox-1_C</td>\n      <td>PF12414.3</td>\n      <td>93</td>\n      <td>sp|O43251|RFOX2_HUMAN</td>\n      <td>-</td>\n      <td>390</td>\n      <td>3.200000e-39</td>\n      <td>133.2</td>\n      <td>29.5</td>\n      <td>1</td>\n      <td>2</td>\n      <td>2.300000e-01</td>\n      <td>6.700000e+02</td>\n      <td>0.7</td>\n      <td>0.0</td>\n      <td>14</td>\n      <td>48</td>\n      <td>177</td>\n      <td>213</td>\n      <td>166</td>\n      <td>243</td>\n      <td>0.66</td>\n      <td>Calcitonin gene-related peptide regulator C te...</td>\n    </tr>\n    <tr>\n      <th>1</th>\n      <td>Fox-1_C</td>\n      <td>PF12414.3</td>\n      <td>93</td>\n      <td>sp|O43251|RFOX2_HUMAN</td>\n      <td>-</td>\n      <td>390</td>\n      <td>3.200000e-39</td>\n      <td>133.2</td>\n      <td>29.5</td>\n      <td>2</td>\n      <td>2</td>\n      <td>8.900000e-42</td>\n      <td>2.600000e-38</td>\n      <td>130.2</td>\n      <td>27.3</td>\n      <td>2</td>\n      <td>93</td>\n      <td>265</td>\n      <td>362</td>\n      <td>264</td>\n      <td>362</td>\n      <td>0.97</td>\n      <td>Calcitonin gene-related peptide regulator C te...</td>\n    </tr>\n    <tr>\n      <th>2</th>\n      <td>RRM_1</td>\n      <td>PF00076.17</td>\n      <td>70</td>\n      <td>sp|O43251|RFOX2_HUMAN</td>\n      <td>-</td>\n      <td>390</td>\n      <td>8.000000e-19</td>\n      <td>67.0</td>\n      <td>0.1</td>\n      <td>1</td>\n      <td>1</td>\n      <td>5.900000e-22</td>\n      <td>1.700000e-18</td>\n      <td>65.9</td>\n      <td>0.1</td>\n      <td>2</td>\n      <td>70</td>\n      <td>124</td>\n      <td>191</td>\n      <td>123</td>\n      <td>191</td>\n      <td>0.97</td>\n      <td>RNA recognition motif. (a.k.a. RRM, RBD, or RN...</td>\n    </tr>\n    <tr>\n      <th>3</th>\n      <td>RRM_6</td>\n      <td>PF14259.1</td>\n      <td>70</td>\n      <td>sp|O43251|RFOX2_HUMAN</td>\n      <td>-</td>\n      <td>390</td>\n      <td>2.400000e-15</td>\n      <td>56.2</td>\n      <td>0.1</td>\n      <td>1</td>\n      <td>1</td>\n      <td>1.400000e-18</td>\n      <td>4.300000e-15</td>\n      <td>55.4</td>\n      <td>0.1</td>\n      <td>1</td>\n      <td>70</td>\n      <td>123</td>\n      <td>191</td>\n      <td>123</td>\n      <td>191</td>\n      <td>0.95</td>\n      <td>RNA recognition motif (a.k.a. RRM, RBD, or RNP...</td>\n    </tr>\n    <tr>\n      <th>4</th>\n      <td>RRM_5</td>\n      <td>PF13893.1</td>\n      <td>56</td>\n      <td>sp|O43251|RFOX2_HUMAN</td>\n      <td>-</td>\n      <td>390</td>\n      <td>8.100000e-11</td>\n      <td>41.6</td>\n      <td>0.1</td>\n      <td>1</td>\n      <td>1</td>\n      <td>5.900000e-14</td>\n      <td>1.800000e-10</td>\n      <td>40.5</td>\n      <td>0.1</td>\n      <td>1</td>\n      <td>54</td>\n      <td>137</td>\n      <td>193</td>\n      <td>137</td>\n      <td>195</td>\n      <td>0.90</td>\n      <td>RNA recognition motif. (a.k.a. RRM, RBD, or RN...</td>\n    </tr>\n  </tbody>\n</table>


### Splice site scoring

If you want to use [MaxEntScan](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) to score splice sites, this is  the feature for you!

To use this feature, you must have [`BioPerl` installed](http://bioperl.org/wiki/Installing_BioPerl).

Say you want to test the strength of these 5' splice sites, which we'll call
`dummy5.fasta`. Notice these 5' splice site sequences are exactly 9 bases long,
as described in the [MaxEntScan `score5ss` documentation](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html).

```
> dummy1
cagGTAAGT
> dummy2
gagGTAAGT
> dummy3
taaATAAGT
```

Then we can get the splice site strength using `score_splice_fasta`

```
from poshsplice.splicestrength import score_splice_fasta
scores = score_splice_fasta('dummy5.fasta', splice_site=5)
```

Notice the `splice_site=5` keyword argument. If you are scoring 3' splice sites
(with length 23 sequences as described in the [MaxEntScan `score3ss` documentation](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html),
then set `splice_site=3`.

The output `scores` is the exact `stdout` output from the original MaxEntScan perl
program, which is a string that looks like this:

```
cagGTAAGT   10.86
gagGTAAGT   11.08
taaATAAGT   0.12
```