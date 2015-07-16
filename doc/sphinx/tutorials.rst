Tutorials
=========

Simple distance estimation
--------------------------

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

Clustering RefSeq
-----------------

Since Refseq Complete contains assemblies and multi-chromosomal organisms that
are not separated by genome, the sequences must first be collated using these
scripts:

:download:`createGiTable.pl <createGiTable.pl>`

:download:`collateRefseqComplete.pl <createGiTable.pl>`

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
:code:`mash sketch` as command line arguments, and create a combined sketch
from the list using list input (:code:`-l`), a k-mer size of 16 (:code:`-k 16`),
and a sketch size of 400 (:code:`-s 400`), which corresponds to an error bound
of 0.05 (this step can be parallelized by dividing the list and combining sketch
files with :code:`mash paste`):

.. code::

  ls refseq*.fna > ../refseq.ls
  mash sketch -l -o ../refseq -k 16 -s 400 ../refseq.ls

From the original directory, run :code:`mash dist` with the combined sketch file
as the reference and query to estimate pairwise distances, using log scaling
(:code:`-L`), 16 threads (:code:`-p 16`), a maximum distance before log-saling
of 0.9 (:code:`-d 0.9`) and a maximum p-value of 1e-10 (:code:`-v 0.0000000001`) for
processing with a clustering tool (and to keep output size manageable):

.. code::

  cd ../
  mash dist -L -d 0.9 -v 0.0000000001 -t -p 16 refseq.msh refseq.msh > distances.tab
  
Cluster base on the log-scaled, pairwise distance estimates, in this case using
`Spici <http://compbio.cs.princeton.edu/spici/>`_:

.. code::

  spici -i distances.tab -o clusters.txt -m 2
