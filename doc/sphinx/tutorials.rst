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

Download the pre-sketched RefSeq archive (reads not provided here;
10x-100x coverage of a single genome with any sequencing technology should
work):

.. download::

`refseq.genomes.k21.s1000.msh <https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh>`_

Plasmids also available:

.. download::

`refseq.genomes+plasmids.k21.s1000.msh <https://gembox.cbcb.umd.edu/mash/refseq.genomes%2Bplasmid.k21s1000.msh>`_
`refseq.plasmids.k21.s1000.msh <https://gembox.cbcb.umd.edu/mash/refseq.plasmid.k21s1000.msh>`_

Concatenate paired ends (this could also be piped to :code:`mash` to save space by
specifying :code:`-` for standard input, zipped or unzipped):

.. code::

 cat reads_1.fastq reads_2.fastq > reads.fastq
 
Sketch the reads, using :code:`-m 2` to improve results
by ignoring single-copy k-mers, which are more likely to be erroneous:

.. code::

  mash sketch -m 2 reads.fastq

Run :code:`mash dist` with the RefSeq archive as the reference and the read
sketch as the query:

.. code::

  mash dist refseq.genomes.k21.s1000.msh reads.fastq.msh > distances.tab

Sort the results to see the top hits and their p-values:

.. code ::

  sort -gk3 distances.tab | head

Screening a read set for containment of RefSeq genomes
------------------------------------------------------

(new in `Mash v2.0 <https://github.com/marbl/Mash/releases>`_)

If a read set potentially has multiple genomes, it can be "screened" against the
database to estimate how well each genome is contained in the read set. We can
use the `SRA Toolkit <https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/>`_ to
download ERR024951:

.. code::

  fastq-dump ERR024951

...and screen it against Refseq Genomes (link above), sorting the results:

.. code::

  mash screen refseq.genomes.k21s1000.msh ERR024951.fastq > screen.tab
  sort -gr screen.tab | head

We see the expected organism, *Salmonella enterica*, but also an apparent contaminant, *Klebsiella pneumoniae*. The fields are [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]:

.. code::

  0.99957	991/1000	26	0	GCF_000841985.1_ViralProj14228_genomic.fna.gz	NC_004313.1 Salmonella phage ST64B, complete genome
  0.99957	991/1000	24	0	GCF_002054545.1_ASM205454v1_genomic.fna.gz	[57 seqs] NZ_MYON01000010.1 Salmonella enterica strain BCW_4905 NODE_10_length_152932_cov_1.77994, whole genome shotgun sequence [...]
  0.999522	990/1000	102	0	GCF_900086185.1_12082_4_85_genomic.fna.gz	[51 seqs] NZ_FLIP01000001.1 Klebsiella pneumoniae strain k1037, whole genome shotgun sequence [...]
  0.999329	986/1000	24	0	GCF_002055205.1_ASM205520v1_genomic.fna.gz	[72 seqs] NZ_MYOO01000010.1 Salmonella enterica strain BCW_4904 NODE_10_length_177558_cov_3.07217, whole genome shotgun sequence [...]
  0.999329	986/1000	24	0	GCF_002054075.1_ASM205407v1_genomic.fna.gz	[88 seqs] NZ_MYNK01000010.1 Salmonella enterica strain BCW_4936 NODE_10_length_177385_cov_3.78874, whole genome shotgun sequence [...]
  0.999329	986/1000	24	0	GCF_000474475.1_CFSAN001184_01.0_genomic.fna.gz	[45 seqs] NZ_AUQM01000001.1 Salmonella enterica subsp. enterica serovar Typhimurium str. CDC_2009K1158 isolate 2009K-1158 SEET1158_1, whole genome shotgun sequence [...]
  0.999329	986/1000	24	0	GCF_000474355.1_CFSAN001186_01.0_genomic.fna.gz	[46 seqs] NZ_AUQN01000001.1 Salmonella enterica subsp. enterica serovar Typhimurium str. CDC_2009K1283 isolate 2009K1283 (Typo) SEET1283_1, whole genome shotgun sequence [...]
  0.999329	986/1000	24	0	GCF_000213635.1_ASM21363v1_genomic.fna.gz	[2 seqs] NC_016863.1 Salmonella enterica subsp. enterica serovar Typhimurium str. UK-1, complete genome [...]
  0.999281	985/1000	24	0	GCF_001271965.1_Salmonella_enterica_CVM_N43825_v1.0_genomic.fna.gz	[67 seqs] NZ_LIMN01000001.1 Salmonella enterica subsp. enterica serovar Typhimurium strain CVM N43825 N43825_contig_1, whole genome shotgun sequence [...]
  0.999281	985/1000	24	0	GCF_000974215.1_SALF-297-3.id2_v1.0_genomic.fna.gz	[90 seqs] NZ_LAPO01000001.1 Salmonella enterica subsp. enterica serovar Typhimurium strain SALF-297-3 NODE_1, whole genome shotgun sequence [...]

Note, however, that multiple strains of *Salmonella enterica* have good identity. This is because they are each contained well when considered independently. For this reason :code:`mash screen` is not a true classifier. However, we can remove much of the redundancy
for interpreting the results using the winner-take-all strategy (:code:`-w`). And while we're at it, let's throw some more cores at
the task to speed it up (:code:`-p 4`):

.. code::

  mash screen -w -p 4 refseq.genomes.k21s1000.msh ERR024951.fastq > screen.tab
  sort -gr screen.tab | head

The output is now much cleaner, with just the two whole genomes, plus phages (a lot of other hits to viruses and assembly contigs would appear further down):

.. code::

  0.99957	991/1000	24	0	GCF_002054545.1_ASM205454v1_genomic.fna.gz	[57 seqs] NZ_MYON01000010.1 Salmonella enterica strain BCW_4905 NODE_10_length_152932_cov_1.77994, whole genome shotgun sequence [...]
  0.99899	979/1000	26	0	GCF_000841985.1_ViralProj14228_genomic.fna.gz	NC_004313.1 Salmonella phage ST64B, complete genome
  0.998844	976/1000	101	0	GCF_900086185.1_12082_4_85_genomic.fna.gz	[51 seqs] NZ_FLIP01000001.1 Klebsiella pneumoniae strain k1037, whole genome shotgun sequence [...]
  0.923964	190/1000	40	0	GCF_000900935.1_ViralProj181984_genomic.fna.gz	NC_019545.1 Salmonella phage SPN3UB, complete genome
  0.900615	111/1000	100	0	GCF_001876675.1_ASM187667v1_genomic.fna.gz	[137 seqs] NZ_MOXK01000132.1 Klebsiella pneumoniae strain AWD5 Contig_(1-18003), whole genome shotgun sequence [...]
  0.887722	82/1000	31	3.16322e-233	GCF_001470135.1_ViralProj306294_genomic.fna.gz	NC_028699.1 Salmonella phage SEN34, complete genome
  0.873204	58/1000	22	1.8212e-156	GCF_000913735.1_ViralProj227000_genomic.fna.gz	NC_022749.1 Shigella phage SfIV, complete genome
  0.868675	52/1000	57	6.26251e-138	GCF_001744215.1_ViralProj344312_genomic.fna.gz	NC_031129.1 Salmonella phage SJ46, complete genome
  0.862715	45/1000	1	1.05185e-116	GCF_001882095.1_ViralProj353688_genomic.fna.gz	NC_031940.1 Salmonella phage 118970_sal3, complete genome
  0.856856	39/1000	21	6.70643e-99	GCF_000841165.1_ViralProj14230_genomic.fna.gz	NC_004348.1 Enterobacteria phage ST64T, complete genome
