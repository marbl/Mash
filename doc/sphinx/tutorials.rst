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

Download and gunzip the pre-sketched RefSeq archive (reads not provided here;
10x-100x coverage of a single genome with any sequencing technology should
work):

.. download::

`RefSeqSketches.msh.gz <http://gembox.cbcb.umd.edu/mash/RefSeqSketches.msh.gz>`_

Concatenate paired ends (this could also be piped to :code:`mash` to save space by
specifying :code:`-` for standard input, zipped or unzipped):

.. code::

 cat reads_1.fastq read_2.fastq > reads.fastq
 
Sketch the reads, using :code:`-m 2` to improve results
by ignoring single-copy k-mers, which are more likely to be erroneous:

.. code::

  mash sketch -m 2 -k 16 -s 400 reads.fastq

Run :code:`mash dist` with the RefSeq archive as the reference and the read
sketch as the query:

.. code::

  mash dist RefSeqSketches.msh reads.fastq.msh > distances.tab

Sort the results to see the top hits and their p-values:

.. code ::

  sort -gk3 distances.tab | head

Containment of RefSeq genomes in a read set
-------------------------------------------

If a read set potentially has multiple genomes, it can be "screened" against the
database to estimate how well each genome is contained in the read set:

.. code::

  mash screen RefSeqSketches.msh reads.fastq > screen.tab

Screening ERR024951, for example, shows the expected organism, *Salmonella enterica*, but also an apparent contaminant, *Klebsiella pneumoniae*, and many contained plasmids (output sorted and abridged; the fields are [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]):

.. code::

  1	400/400	93	0	./rcn/refseq-NG-623-.-.-.-pSTR1-Shigella_flexneri.fna	
  1	400/400	78	0	./rcn/refseq-NZ-573-PRJNA224116-SAMN02885365-GCF_000764615.1-pPMK1_A-Klebsiella_pneumoniae.fna	
  1	400/400	76	0	./rcn/refseq-NG-1351-.-.-.-pNCC801-Enterococcus_faecalis.fna	
  1	400/400	54	0	./rcn/refseq-NG-98360-.-.-.-virulence_plasmid-Salmonella_enterica_subsp._enterica_serovar_Dublin.fna
  [...]
  0.999844	399/400	1130	0	./rcn/refseq-NR-90371-PRJNA188943-.-.-.-Salmonella_enterica_subsp._enterica_serovar_Typhimurium.fna
  [...]
  0.998897	393/400	113	0	./rcn/refseq-NZ-1328378-PRJNA224116-SAMN02138586-GCF_000492775.1-.-Klebsiella_pneumoniae_MGH_46.fna
  [...]

Building a custom RefSeq database
---------------------------------

To create the RefSeq Mash database, genomes were downloaded from NCBI
(:code:`ftp.ncbi.nlm.nih.gov/refseq/release/complete`, fasta sequence and
GenBank annotations for :code:`genomic`), and the
`refseqCollate <https://github.com/ondovb/refseqCollate/releases>`_ utility was
used to collate contigs/chromosomes into individual fasta files per genome.
Groups of these files were sketched in parallel and then pasted together with
:code:`mash paste`. This process could be repeated for more current or custom
databases.
