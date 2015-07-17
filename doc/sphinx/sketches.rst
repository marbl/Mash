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
compared with :code:`mash dist`. A good rule for k-mer size is:

.. math::
  \log_4 {G}
  
...where :math:`G` is the average size of genomes that will be compared.

sketch size
'''''''''''
Sketch size corresponds to the number of (non-redundant) min-hashes that are
kept. Larger sketches will better represent the sequence, but at the cost of
larger sketch files and longer comparison times. The error bound of a distance
estimation for a given sketch size :math:`s` is formulated as:

.. math::
  \sqrt(\frac{1}{s})

Sketch size is specified with :code:`-s`. Sketches of different sizes can be
compared with :code:`mash dist`, although the comparison will be restricted to
the smaller of the two sizes.

Strand-independence with canonical k-mers
-----------------------------------------
By default, :code:`mash` will ignore strandedness when sketching by using
canonical k-mers, as done in `Jellyfish`_. This works by using the reverse
complement of a k-mer if it comes before the original k-mer alphabetically.
It also means k-mers that do not contain only nucleotides (A, C, G, T, and their
lowercases) must be ignored. To use every k-mer as it appears, :code:`-n`
(noncanonical) can be specified when sketching.

Cleaning up read sets with Bloom filtering
------------------------------------------

Since MinHash is a k-mer based method, removing unique k-mers greatly improves
results for read sets, since unique k-mers are likely to represent sequencing
error. :code:`mash` provides an efficient way to filter without prior k-mer
counting by using a Bloom filter. This method can underfilter, but it will
never overfilter (non-unique k-mers are guaranteed to be kept), and it requires
significantly less time and memory than true k-mer counting. The filter can be
enabled with :code:`-u` when sketching (in :code:`mash sketch` or :code:`mash
dist`). The amount of underfiltering can be managed with the parameters of the
Bloom filter (:code:`-g`, :code:`-e`, and :code:`-m`).

Working with sketch files
-------------------------

The sketch or sketches stored in a sketch file, and their parameters, can be 
inspected with :code:`mash info`. If sketch files have matching k-mer sizes,
their sketches can be combined into a single file with :code:`mash paste`. This
allows simple pairwise comparisons with :code:`mash dist`, and allows sketching
of multiple files to be parallelized.

.. _Jellyfish: http://www.cbcb.umd.edu/software/jellyfish/
