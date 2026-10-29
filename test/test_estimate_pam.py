import random
import types
import unittest

import numpy as np
import pyopa

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def random_alignment(rng, length, subst_rate, gap_rate):
    """returns two aligned rows of equal length with random substitutions and gaps"""
    a = ''.join(rng.choice(AMINO_ACIDS) for _ in range(length))
    b = ''.join('_' if rng.random() < gap_rate else
                (rng.choice(AMINO_ACIDS) if rng.random() < subst_rate else c) for c in a)
    a = ''.join('_' if rng.random() < gap_rate / 2 else c for c in a)
    return pyopa.Sequence(a), pyopa.Sequence(b)


class EstimatePamTest(unittest.TestCase):
    """estimate_pam (eigendecomposition of logPAM1) against the original EstimatePam"""

    def setUp(self):
        self.defaults = pyopa.load_default_environments()
        self.dms = pyopa.MutipleAlEnv(self.defaults['environments'], self.defaults['log_pam1'])

    def assertMatchesReference(self, dms, s1, s2):
        res = dms.estimate_pam(s1, s2)
        ref = dms._estimate_pam_reference(s1, s2)
        msg = '%s / %s: %s != %s' % (s1.convert_readable(), s2.convert_readable(), res, ref)
        # the score and the pam distance come from the same search over the DMS
        self.assertEqual(res[:2], ref[:2], msg)
        # the variance is computed with a different, but equivalent, formula
        np.testing.assert_allclose(res[2], ref[2], rtol=1e-8, err_msg=msg)

    def test_random_alignments(self):
        rng = random.Random(7)
        for _ in range(300):
            length = rng.choice([10, 50, 200, 1000])
            s1, s2 = random_alignment(rng, length, subst_rate=rng.random() ** 2,
                                      gap_rate=rng.random() * 0.1)
            self.assertMatchesReference(self.dms, s1, s2)

    def test_edge_cases(self):
        cases = [
            ('ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVWY'),  # identical, smallest pam
            ('A', 'A'),
            ('A', 'C'),
            ('AAAAAAAAAA', 'CCCCCCCCCC'),  # saturated, largest pam
            ('ACDE__GHIK', 'AC__FGHIKL'),
            ('ACDEXXGHIK', 'ACDEFGXXIK'),  # X is not counted
            ('ABZDEFGHIK', 'ABZDEFGHIK'),  # ambiguous symbols
            ('___', 'ACD'),
            ('A_C', 'A_C'),
        ]
        for a, b in cases:
            self.assertMatchesReference(self.dms, pyopa.Sequence(a), pyopa.Sequence(b))

    def test_estimate_many_pam(self):
        rng = random.Random(3)
        ancestor = ''.join(rng.choice(AMINO_ACIDS) for _ in range(300))
        seqs = [pyopa.Sequence(ancestor)] + [pyopa.Sequence(''.join(
            '_' if rng.random() < 0.02 else (rng.choice(AMINO_ACIDS) if rng.random() < 0.1 * i else c)
            for c in ancestor)) for i in range(1, 8)]
        score, dist, var = self.dms.estimate_many_pam(seqs, num_threads=2)
        k = 0
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                self.assertEqual([score[k], dist[k], var[k]], self.dms.estimate_pam(seqs[i], seqs[j]))
                k += 1

    def test_non_reversible_log_pam1_falls_back(self):
        # breaking detailed balance of logPAM1 (columns still sum to 0)
        # makes estimate_pam use the original EstimatePam
        q = self.defaults['log_pam1'].float64_matrix.copy()
        a, b = (ord(c) - ord('A') for c in 'AC')
        q[b, b] -= 0.5 * q[a, b]
        q[a, b] *= 1.5
        # MutipleAlEnv only uses the float64_matrix of log_pam1
        dms = pyopa.MutipleAlEnv(self.defaults['environments'], types.SimpleNamespace(float64_matrix=q))

        rng = random.Random(5)
        seqs = [random_alignment(rng, 200, subst_rate=0.3, gap_rate=0.05) for _ in range(20)]
        for s1, s2 in seqs:
            self.assertEqual(dms.estimate_pam(s1, s2), dms._estimate_pam_reference(s1, s2))
        score, dist, var = dms.estimate_many_pam(list(seqs[0]))
        self.assertEqual([score[0], dist[0], var[0]], dms._estimate_pam_reference(*seqs[0]))


if __name__ == '__main__':
    unittest.main()
