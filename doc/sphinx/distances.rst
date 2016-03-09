Distance Estimation
===================

MinHash Jaccard estimation
--------------------------

Given :math:`k`-mer sets :math:`A` and :math:`B`, the MinHash algorithm provides an
estimation of the Jaccard index:

.. math::

 j(A_s,B_s) = \frac {\lvert A_s \cap B_s \rvert} s

where :math:`A_s` and :math:`B_s` are subsets such that
:math:`\lvert A_s \cup B_s \rvert` is equal to the sketch size, :math:`s`,
allowing for a known error bound as suggested by Broder [#f1]_. This is done by
using a merge-sort algorithm to find common values between the two sorted
sketches and terminating when the total number of hashes seen reaches the sketch
size (or all hashes in both sketches have been seen).

Mash distance formulation
-------------------------

For mutating a sequence with :math:`t` total :math:`k`-mers and a
conserved :math:`k`-mer count :math:`w`, an approximate mutation rate :math:`d` can be
estimated using a Poisson model of mutations occurring in :math:`k`-mers, as suggested
by Fan et al. [#f2]_:

.. math::

  d = \frac {-1} k \ln \frac w t

In order to use a Jaccard estimate :math:`j` between two :math:`k`-mer sets of arbitrary
sizes, the Jaccard estimate can be framed in terms of the conserved :math:`k`-mer count
:math:`w` and the average set size :math:`n`:

.. math::

  j \approx \frac w {2n - w}

To substitute :math:`n` for the total :math:`k`-mer count :math:`t` in the mutation
estimation, this approximation can be reformulated as:

.. math::

  \frac w n \approx \frac {2j} {1 + j}

Substituting :math:`\frac w n` for :math:`\frac w t` thus yields the Mash
distance:

.. math::

  D(k,j)=\Bigg\{\begin{split}
  &1&,\ j=0\\
  &\frac {-1} k \ln \frac {2j} {1 + j}&,\ 0<j\le 1
  \end{split}



  
Assessing significance with p-values
------------------------------------
Since MinHash distances are probabilistic estimates, it is important to
consider the probability of seeing a given distance by chance. :code:`mash dist`
thus provides p-values with distance estimations. Lower p-values correspond to
more confident distance estimations, and will often be rounded down to 0 due to
floating point limits. If p-values are high (above, say, 0.01), the :math:`k`-mer size
is probably too small for the size of the genomes being compared.

When estimating the distance of genome 1 and genome 2 from sketches with the
properties:

  :math:`\Sigma` := alphabet
  
  :math:`k` := :math:`k`-mer size
  
  :math:`l_1` := length of genome 1
  
  :math:`l_2` := length of genome 2
  
  :math:`s` := sketch size
  
  :math:`x` := number of shared :math:`k`-mers between sketches of size :math:`s` of
  genome 1 and genome 2
  
...the chance of a :math:`k`-mer appearing in random sequences of lengths :math:`l_1`
and :math:`l_2` are estimated as:

.. math::

  r_1 = 1-(1-\frac{1}{{\lvert\Sigma\rvert}^k})^{l_1} \approx \frac{l_1}{l_1+{\lvert\Sigma\rvert}^k}
  
  r_2 = 1-(1-\frac{1}{{\lvert\Sigma\rvert}^k})^{l_2} \approx \frac{l_2}{l_2+{\lvert\Sigma\rvert}^k}
  
The expected Jaccard index of the sketches of the random sequences is then:

.. math::

  j_r = \frac{r_1 r_2}{r_1 + r_2 - r_1 r_2}

...and the probability of observing at least :math:`x` shared :math:`k`-mers can be
estimated with the tail of a cumulative binomial distribution:

.. math::
  
  p = 1 - \sum\limits_{i=0}^{x-1} \binom{s}{i} j_r^i (1-j_r)^{s-i}

.. [#f1] `Broder, A.Z. On the resemblance and containment of documents. Compression and Complexity of Sequences 1997 - Proceedings, 21-29 (1998). <http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=666900&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D666900>`_
.. [#f2] `Fan, H., Ives, A.R., Surget-Groba, Y. & Cannon, C.H. An assembly and alignment-free method of phylogeny reconstruction from next-generation sequencing data. BMC genomics 16, 522 (2015). <http://www.ncbi.nlm.nih.gov/pubmed/26169061>`_