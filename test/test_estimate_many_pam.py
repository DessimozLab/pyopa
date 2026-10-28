import itertools
import os
import random
import unittest

import numpy as np
import pyopa


def mutate(seq, rng, subst_rate, gap_rate, alphabet='ACDEFGHIKLMNPQRSTVWY'):
    """returns a copy of seq with random substitutions and gaps, keeping its length"""
    res = []
    for c in seq:
        r = rng.random()
        if r < gap_rate:
            res.append('_')
        elif r < gap_rate + subst_rate:
            res.append(rng.choice(alphabet))
        else:
            res.append(c)
    return ''.join(res)


class EstimateManyPamTest(unittest.TestCase):

    def setUp(self):
        defaults = pyopa.load_default_environments()
        self.dms = pyopa.MutipleAlEnv(defaults['environments'], defaults['log_pam1'])

        with open(os.path.join(os.path.dirname(__file__), 'data', 'testseqs.txt')) as f:
            ancestor = f.readline().strip()[:400]

        # rows of an artificial multiple sequence alignment, all of equal length
        # and with increasing divergence from the ancestor
        rng = random.Random(42)
        rows = [ancestor] + [mutate(ancestor, rng, subst_rate=0.05 * i, gap_rate=0.02)
                             for i in range(1, 7)]
        self.seqs = [pyopa.Sequence(r) for r in rows]

    def expected(self, seqs):
        """reference results with estimate_pam, in condensed (upper triangle) pair order"""
        res = [self.dms.estimate_pam(s1, s2) for s1, s2 in itertools.combinations(seqs, 2)]
        return tuple(np.array(col) for col in zip(*res))

    def test_matches_estimate_pam(self):
        score, dist, var = self.dms.estimate_many_pam(self.seqs)
        exp_score, exp_dist, exp_var = self.expected(self.seqs)

        self.assertEqual(len(score), len(self.seqs) * (len(self.seqs) - 1) // 2)
        np.testing.assert_allclose(score, exp_score, rtol=1e-10, err_msg='similarity scores differ')
        np.testing.assert_allclose(dist, exp_dist, rtol=1e-10, err_msg='pam distances differ')
        np.testing.assert_allclose(var, exp_var, rtol=1e-10, err_msg='variances differ')

    def test_two_sequences(self):
        score, dist, var = self.dms.estimate_many_pam(self.seqs[:2])
        ref = self.dms.estimate_pam(self.seqs[0], self.seqs[1])
        np.testing.assert_allclose([score[0], dist[0], var[0]], ref, rtol=1e-10)

    def test_independent_of_num_threads(self):
        # many short sequences give many quick estimates, which makes
        # races between threads likely to show up within a few repetitions
        rng = random.Random(1)
        ancestor = self.seqs[0].convert_readable()[:20]
        seqs = [pyopa.Sequence(mutate(ancestor, rng, subst_rate=0.3, gap_rate=0.02)) for _ in range(150)]
        single = self.dms.estimate_many_pam(seqs, num_threads=1)
        for _ in range(10):
            multi = self.dms.estimate_many_pam(seqs, num_threads=8)
            for a, b in zip(single, multi):
                np.testing.assert_array_equal(a, b)

    def test_invalid_input(self):
        self.assertRaises(ValueError, self.dms.estimate_many_pam, [])
        self.assertRaises(ValueError, self.dms.estimate_many_pam, self.seqs[:1])
        shorter = pyopa.Sequence(self.seqs[1].convert_readable()[:-1])
        self.assertRaises(ValueError, self.dms.estimate_many_pam, [self.seqs[0], shorter])


if __name__ == '__main__':
    unittest.main()
