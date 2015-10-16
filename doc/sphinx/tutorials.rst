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

