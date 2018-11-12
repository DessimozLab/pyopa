import gzip
import json
import time
import pyopa
import os
import unittest
# google-perfprofiler wrapper
#import yep


class RunTimeOfLongSequenceTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with gzip.open(
                os.path.join(os.path.dirname(__file__),
                             "data", "long_seq_pair.json.gz"),
                'rt') as fh:
            cls.data = json.load(fh)
        cls.aligner = pyopa.Aligner()

    def test_runtime_single_matrix(self):
        env = self.aligner.environment_at_distance(self.data['distance'])
        s1 = pyopa.Sequence(str(self.data['s1']))
        s2 = pyopa.Sequence(str(self.data['s2']))
        t0 = time.time()
        nr_runs = 5
        for x in range(nr_runs):
            #yep.start('align_{}.prof'.format(x))
            double_alignment = pyopa.align_double(s1, s2, env, stop_at_threshold=False,
                                                  is_global=False, calculate_ranges=True)
            as1, as2 = pyopa.align_strings(s1, s2, env, False, double_alignment)
            #yep.stop()
        print("Avg time used to compute alignment on fixed matrix: {}sec"
              .format((time.time()-t0) / nr_runs))
        print("Darwin's run time for this alignment: {}sec".format(self.data['time_single_matrix_align']))
        self.assertEqual(str(self.data['as1']), as1)
        self.assertEqual(str(self.data['as2']), as2)

        #self.aligner.align(self.data['s1'], self.data['s2'], env=env)


if __name__ == '__main__':
    unittest.main()
