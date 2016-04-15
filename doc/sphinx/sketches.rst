Sketches
========

For sequences to be compared with :code:`mash`, they must first be `sketched`,
which creates vastly reduced representations of them. This will happen
automatically if :code:`mash dist` is given raw sequences. However, if multiple
comparisons will be performed, it is more efficient to create sketches with
:code:`mash sketch` first and provide them to :code:`mash dist` in place of the
raw sequences. Sketching parameters can be provided to either tool via
command line options.

Reduced representations with MinHash tables
-------------------------------------------
Sketches are used by the `MinHash` algorithm to allow fast distance estimations
with low storage and memory requirements. To make a sketch, each k-mer in a
sequence is `hashed`, which creates a pseudo-random identifier. By sorting these
identifiers (`hashes`), a small subset from the top of the sorted list can
represent the entire sequence (these are `min-hashes`). The more similar another
sequence is, the more min-hashes it is likely to share.

k-mer size
''''''''''
As in any k-mer based method, larger k-mers will provide more specificity, while
smaller k-mers will provide more sensitivity. Larger genomes will also require
larger k-mers to avoid k-mers that are shared by chance. K-mer size is
specified with :code:`-k`, and sketch files must have the same k-mer size to be
compared with :code:`mash dist`. When :code:`mash sketch` is run, it
automatically assesses the specified k-mer size against the sizes of input
genomes by estimating the probability of a random match as:

.. math::
  p = \frac 1 {\frac {\left(\overline\Sigma\right)^k} g + 1}
  
...where :math:`g` is the genome size and :math:`\Sigma` is the alphabet (ACGT
by default). If this probability exceeds a threshold (specified by
:code:`-w`; 0.01 by default) for any input genomes, a warning will be given
with the minimum k-mer size needed to get within the threshold.

For large collections of sketches, memory and storage may also be a
consideration when choosing a k-mer size. Mash will use 32-bit hashes, rather
than 64-bit, if they can encompass the full k-mer space for the alphabet in use.
This will (roughly) halve the size of the size of the sketch file on disk and
the memory it uses when loaded for :code:`mash dist`. The criterion for using a
32-bit hash is:

.. math::
   \left({\overline\Sigma}\right)^k \leq 2^{32}

...which becomes :math:`k \leq 16` for nucleotides (the default) and
:math:`k \leq 7` for amino acids.

sketch size
'''''''''''
Sketch size corresponds to the number of (non-redundant) min-hashes that are
kept. Larger sketches will better represent the sequence, but at the cost of
larger sketch files and longer comparison times. The error bound of a distance
estimation for a given sketch size :math:`s` is formulated as:

.. math::
  \sqrt{\frac{1}{s}}

Sketch size is specified with :code:`-s`. Sketches of different sizes can be
compared with :code:`mash dist`, although the comparison will be restricted to
the smaller of the two sizes.

Strand and alphabet
-------------------
By default, :code:`mash` uses a nucleotide alphabet (ACGT), is case-insensitive,
and will ignore strandedness by using canonical k-mers, as done in
`Jellyfish`_. This works by using the reverse complement of a k-mer if it comes
before the original k-mer alphabetically. Strandedness can be preserved with
:code:`-n` (non-canonical) and case can be preserved with :code:`-Z`. Note that
the default nucleotide alphabet does not include lowercase and thus will filter
out k-mers with lowercase nucleotides if :code:`-Z` is specified. The amino acid
alphabet can be specified with :code:`-a`, which also changes the default k-mer
size to reflect the denser information. A completely custom alphabet can also be
specified with :code:`-z`. Note that alphabet size affects p-value calculation
and hash size (see `Assessing significance with p-values <distances.html#assessing-significance-with-p-values>`_ and `k-mer size`_).


Sketching read sets
-------------------

When sketching reads instead of complete genomes or assemblies, :code:`-r`
should be specified, which will estimate genome size from k-mer content
rather than total sequence length, allowing more accurate p-vlaues. Genome
size can also be specified directly with :code:`-g`. Additionally, Since
MinHash is a k-mer based method, removing unique or low-copy k-mers usually
improves results for read sets, since these k-mers are likely to represent
sequencing error. The minimum copies of each k-mer required can be specified
with :code:`-m` (e.g. :code:`-m 2` to filter unique). However, this could
lead to high memory usage if genome size is high and coverage is low, such as
in metagenomic read sets. In these cases a Bloom filter can be used (:code:`-b`)
to filter out most unique k-mers with constant memory. If coverage is high (e.g.
>100x), it can be helpful to limit it to save time and to avoid repeat errors
appearing as legitimate k-mers. This can be done with :code:`-c`, which stops
sketching reads once the estimated average coverage (based on k-mer
multiplicity) reaches the target.

Working with sketch files
-------------------------

The sketch or sketches stored in a sketch file, and their parameters, can be 
inspected with :code:`mash info`. If sketch files have matching k-mer sizes,
their sketches can be combined into a single file with :code:`mash paste`. This
allows simple pairwise comparisons with :code:`mash dist`, and allows sketching
of multiple files to be parallelized.

.. _Jellyfish: http://www.cbcb.umd.edu/software/jellyfish/
