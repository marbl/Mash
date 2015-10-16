Tutorials
=========

Simple distance estimation
--------------------------

.. download::
`genome1.fna <https://github.com/marbl/Mash/raw/master/data/genome1.fna>`_
`genome2.fna <https://github.com/marbl/Mash/raw/master/data/genome2.fna>`_
`genome3.fna <https://github.com/marbl/Mash/raw/master/data/genome3.fna>`_

.. code::

  mash dist genome1.fna genome2.fna

Saving time by sketching first
------------------------------

.. code::

  mash sketch genome1.fna
  mash sketch genome2.fna
  mash dist genome1.fna.msh genome2.fna.msh

Pairwise comparisons with compound sketch files
-----------------------------------------------

.. code::

  mash sketch -o reference genome1.fna genome2.fna
  mash info reference.msh
  mash dist reference.msh genome3.fna

Querying read sets against an existing RefSeq sketch
----------------------------------------------------

Download the pre-sketched RefSeq archive:

.. download::

`refseq.msh <https://github.com/marbl/Mash/raw/master/data/refseq.msh>`_

Run :code:`mash dist` with the archive as the reference and the read set as the
query, using :code:`-u` to improve results by filtering unique k-mers:

.. code::

  mash dist -u refseq.msh reads.fastq > distances.tab

Sort the results to see the top hits and their p-values:

.. code ::

  sort -nrk3 distances.tab | head

Clustering RefSeq from a fresh download
---------------------------------------

The efficiency of :code:`mash` allows the *de novo* clustering of all genomes in
RefSeq Complete to be performed on a typical workstation in less than a day.

Since Refseq Complete contains assemblies and multi-chromosomal organisms that
are not separated by genome, the sequences must first be collated using these
scripts:

[ Coming soon ]

.. `createGiTable.pl <createGiTable.pl>`

.. `collateRefseqComplete.pl <createGiTable.pl>`

First, create directories for the raw and collated files (there will be tens of
thousands of collated files):

.. code::

  mkdir refseq
  mkdir collated
  
Download the raw fasta and Genbank records into the :code:`refseq` directory:

.. code::

  rsync -a ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.*.genomic.gbff.gz refseq
  rsync -a ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.*.genomic.fna.gz refseq
  
Associate GIs with genomes based on metadata from the Genbank records (writes :code:`giTable.tab`):

.. code::

  ./createGiTable.pl refseq/complete.*.genomic.gbff.gz
  
From the :code:`collated` directory, collate fasta sequences (one genome to a
file) based on the table:

.. code::

  cd collated
  ../collateRefseqComplete.pl ../refseq/complete.*.genomic.fna.gz

Write a list of the collated files, since there are too many to give to
:code:`mash sketch` as command line arguments:

.. code ::

  ls refseq*.fna > ../refseq.ls

From the original directory, create the combined sketch
using list input (:code:`-l`), a k-mer size of 16 (:code:`-k 16`),
and a sketch size of 400 (:code:`-s 400`), which corresponds to an error bound
of 0.05 (this step can be parallelized by dividing the list and combining sketch
files with :code:`mash paste`):

.. code::

  cd ../
  mash sketch -l -o refseq -k 16 -s 400 refseq.ls

Run :code:`mash dist` with the combined sketch file as both the reference and
query to estimate pairwise distances, using 16 threads (:code:`-p 16`), a
maximum distance of 0.05 (:code:`-d 0.05`) and a maximum p-value of 1e-10
(:code:`-v 0.0000000001`).
:

.. code::

  mash dist -d 0.05 -v 0.0000000001 -p 16 refseq.msh refseq.msh > distances.tab
  
Cluster based on the log-scaled, pairwise distance estimates, in this case using
`Spici <http://compbio.cs.princeton.edu/spici/>`_, with :code:`-m 2` to specify
a large, sparse graph to keep memory manageable:

.. code::

  spici -i distances.tab -o clusters.txt -m 2

The output (:code:`clusters.txt`) will have one cluster per line, with genomes
named by :code:`collateRefseqComplete.pl` as follows:

refseq-[PRX]-[TAX]-[PRJ]-[SMP]-[ASM]-[PLS]-[ORG].fna

...where:

- [PRX] = two-letter accession prefix
- [TAX] = taxonomy ID
- [PRJ] = BioProject ID, if available
- [SMP] = BioSample ID, if available
- [ASM] = Assembly ID, if available
- [PLS] = Plasmid ID, if available
- [ORG] = Organism name

(missing fields will be ".")