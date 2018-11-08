PyOPA - Optimal Pairwise Alignments
===================================

This python package provides a fast implementation to compute 

  - optimal pairwise alignments of molecular sequences

  - ML distance estimates of pairwise alignments.

The implementation uses `Farrar's algorithm <http://bioinformatics.oxfordjournals.org/content/23/2/156.abstract>_`
to compute the optimal pairwise alignment using SSE vectorization operations.
This package implements the Smith-Waterman and Needleman-Wunsch algorithm to
compute the local and global sequence alignments.

Example
-------

.. code-block:: python

   import pyopa
   log_pam1_env = pyopa.read_env_json(os.path.join(pyopa.matrix_dir(), 'logPAM1.json'))
   s1 = pyopa.Sequence('GCANLVSRLENNSRLLNRDLIAVKINADVYKDPNAGALRL')
   s2 = pyopa.Sequence('GCANPSTLETNSQLVNRELIAVKINPRVYKGPNLGAFRL')

   # super fast check whether the alignment reaches a given min-score
   min_score = 100
   pam250_env = pyopa.generate_env(log_pam1_env, 250, min_score)
   pyopa.align_short(s1, s2, pam250_env)

   
