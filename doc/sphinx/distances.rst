Distance Estimation
===================

Jaccard distance formulation
----------------------------

For sequences with sketches :math:`A` and :math:`B`, the distance of the
sequences can be estimated by the Jaccard distance of the sketches:

.. math::

 d_J(A,B) = 1 - \frac {\lvert A \cap B \rvert} {\lvert A \cup B \rvert}

As implemented in :code:`mash`, the subsets :math:`A_s` and :math:`B_s` are used
such that :math:`\lvert A_s \cup B_s \rvert` is equal to the sketch size,
allowing for a known error bound based on the given sketch size. This is done by
using a merge-sort algorithm to find common values between the two sorted
sketches and terminating when the total number of hashes seen reaches the sketch
size (or all hashes in both sketches have been seen), as suggested by Broder [#f1]_.

For closer correlation with nucleotide identity, the distance can be log scaled
(:code:`-L` in :code:`mash dist`), using pseudo-counts to avoid taking the log
of 0 when no hashes are shared:

.. math::

  -\log_{10} \frac {\lvert A \cap B \rvert + 1} {\lvert A \cup B \rvert + 1}

Assessing significance with p-values
------------------------------------
Since MinHash distances are probabilistic estimates, it is important to
consider the probability of seeing a given distance by chance. :code:`mash dist`
thus provides p-values with distance estimations. Lower p-values correspond to
more confident distance estimations, and will often be rounded down to 0 due to
floating point limits. If p-values are high (above, say, 0.01), the k-mer size
is probably too small for the size of the genomes being compared.

When estimating the distance of genome 1 and genome 2 from sketches with the
properties:

  :math:`k` := k-mer size
  
  :math:`l_1` := length of genome 1
  
  :math:`l_2` := length of genome 2
  
  :math:`s` := sketch size
  
  :math:`x` := number of shared k-mers between sketches of size :math:`s` of
  genome 1 and genome 2
  
...the chance of a k-mer appearing in random sequences of lengths :math:`l_1`
and :math:`l_2` are estimated as:

.. math::

  r_1 = 1-(1-\frac{1}{4^k})^{l_1} \approx \frac{l_1}{l_1+4^k}
  
  r_2 = 1-(1-\frac{1}{4^k})^{l_2} \approx \frac{l_2}{l_2+4^k}
  
The expected Jaccard index of the sketches of the random sequences is then:

.. math::

  j = \frac{r_1 r_2}{r_1 + r_2 - r_1 r_2}

...and the probability of observing at least :math:`x` shared k-mers can be
estimated with the tail of a cumulative binomial distribution:

.. math::
  
  p = 1 - \sum\limits_{i=0}^{x-1} \binom{s}{i} j^i (1-j)^{s-i}

.. [#f1] Broder ...
