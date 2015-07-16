Distance Estimation
===================

Jaccard distance formulation
----------------------------

Assessing significance with p-values
------------------------------------

  :math:`l_1` := length of genome 1
  
  :math:`l_2` := length of genome 2
  
  :math:`s` := sketch size
  
.. math::

  r_1 = 1-(1-\frac{1}{4^k})^{l_1} \approx \frac{l_1}{l_1+4^k}
  
  r_2 = 1-(1-\frac{1}{4^k})^{l_2} \approx \frac{l_2}{l_2+4^k}
  
  j = \frac{r_1 r_2}{r_1 + r_2 - r_1 r_2}
  
  p = 1 - \sum\limits_{i=0}^{x-1} \binom{s}{i} j^i (1-j)^{s-i}
