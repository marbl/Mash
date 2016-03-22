Tutorials
=========

Simple distance estimation
--------------------------

Download example *E. coli* genomes:

.. download::
| `genome1.fna <https://gembox.cbcb.umd.edu/mash/genome1.fna>`_ 
| `genome2.fna <https://gembox.cbcb.umd.edu/mash/genome2.fna>`_

Run:

.. code::

  mash dist genome1.fna genome2.fna

The results are tab delimited lists of Reference-ID, Query-ID, Mash-distance,
P-value, and Matching-hashes:

.. code::

  genome1.fna	genome2.fna	0.0222766	0	456/1000

Saving time by sketching first
------------------------------

.. code::

  mash sketch genome1.fna
  mash sketch genome2.fna
  mash dist genome1.fna.msh genome2.fna.msh

Pairwise comparisons with compound sketch files
-----------------------------------------------

Download additional example *E. coli* genome:

| `genome3.fna <https://gembox.cbcb.umd.edu/mash/genome3.fna>`_

Sketch the first two genomes to create a combined archive, use :code:`mash info`
to verify its contents, and estimate pairwise distances:

.. code::

  mash sketch -o reference genome1.fna genome2.fna
  mash info reference.msh
  mash dist reference.msh genome3.fna

This will estimate the distance from each query (which there is one of) to each
reference (which there are two of in the sketch file):

.. code::

  genome1.fna	genome3.fna	0	0	1000/1000
  genome2.fna	genome3.fna	0.0222766	0	456/1000

Querying read sets against an existing RefSeq sketch
----------------------------------------------------

Download and gunzip the pre-sketched RefSeq archive:

.. download::

`RefSeqSketches.msh.gz <http://gembox.cbcb.umd.edu/mash/RefSeqSketches.msh.gz>`_

Sketch the reads (not provided here; 10x-100x coverage of a single bacterial genome
with any sequencing technology should work), using :code:`-u` to improve results
by filtering unique k-mers:

.. code::

  mash sketch -u reads.fastq

Run :code:`mash dist` with the RefSeq archive as the reference and the read
sketch as the query:

.. code::

  mash dist RefSeqSketches.msh reads.fastq.msh > distances.tab

Sort the results to see the top hits and their p-values:

.. code ::

  sort -gk3 distances.tab | head
