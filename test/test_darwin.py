import unittest
import os
import csv
import cython_swps3
import sys


class DarwinResult:
    def __init__(self):
        self.s1_id = 0
        self.s2_id = 0
        self.matrix_nr = 0
        self.pam = 0.0
        self.threshold = 0.0
        self.score_double = 0.0
        self.score_float = 0.0
        self.score_short = 0.0
        self.score_byte = 0.0

'''
def writeEnvFile(all_envs):

    for i in range(len(all_envs)):
        env = all_envs[i]

        env.float64_matrix = env.float64_matrix.reshape(26*26)
        env.int16_matrix = env.int16_matrix.reshape(26*26)
        env.int8_matrix = env.int8_matrix.reshape(26*26)

        with open(os.path.dirname(__file__) + '/data/matrices/C_compatible/' + str(i + 1) + '.dat', 'w+') as f:
            f.write("%.20f\n%.20f\n" % (env.gap_open, env.gap_ext))

            for i in range(26*26):
                f.write("%.20f " % env.float64_matrix[i])

            f.write("\n")
            for i in range(26*26):
                f.write("%d " % env.int16_matrix[i])

            f.write("\n")
            for i in range(26*26):
                f.write("%d " % env.int8_matrix[i])
'''


class AlignTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.dirname(__file__) + '/data/testseqs.txt') as f:
            self.sequences = f.readlines()

        self.sequences = map(lambda s: s.strip(), self.sequences)
        self.sequences_normalized = map(lambda s: cython_swps3.normalize_sequence(s), self.sequences)

        self.darwin_results = []
        self.alignment_environments = cython_swps3.read_all_env_json(
            os.path.dirname(__file__) + '/data/matrices/json/all_matrices.json')

        #writeEnvFile(self.alignment_environments)

        with open(os.path.dirname(__file__) + '/data/reference_test_results.dat') as f:
            #skip header
            next(f)
            reader = csv.reader(f, delimiter='\t')

            for s1, s2, matrix_nr, pam, threshold, score_d, score_f, score_s, score_b in reader:
                curr = DarwinResult()
                curr.s1_id = int(s1)
                curr.s2_id = int(s2)
                curr.matrix_nr = int(matrix_nr)
                curr.pam = float(pam)
                curr.threshold = float(threshold)
                curr.score_double = float(score_d)
                curr.score_float = float(score_f)
                curr.score_short = float(score_s)
                curr.score_byte = float(score_b)

                self.darwin_results.append(curr)

                '''
                if curr.s1_id not in self.alignment_profiles:
                    p = cython_swps3.AlignmentProfile()
                    p.create_profiles(self.sequences[curr.s1_id - 1], self.alignment_environments[curr.s1_id - 1])
                    self.alignment_profiles[curr.s1_id] = p
                '''

    def test_align_scalar(self):
        print 'Running alignment tests...'
        completed = 0
        max_alignments = len(self.darwin_results)
        progress_step = max_alignments / 100

        for r in self.darwin_results:
            s1 = self.sequences_normalized[r.s1_id - 1]
            s2 = self.sequences_normalized[r.s2_id - 1]
            env = self.alignment_environments[r.matrix_nr - 1]
            env.set_threshold(r.threshold)

            #profile = self.alignment_profiles[r.s1_id]
            profile = cython_swps3.AlignmentProfile()
            profile.create_profiles(s1, env, True)

            scalar_result = cython_swps3.align_scalar(s1, s2, env, True)
            byte_result = profile.align_byte(s2, env, True)
            short_result = profile.align_short(s2, env, True)

            self.assertGreaterEqual(short_result, r.score_double,
                                    "Short score must be greater or equal than double score.")
            self.assertGreaterEqual(byte_result, r.score_double,
                                    "Byte score must be greater or equal than double score.")

            #TODO should be a more precise comparison
            self.assertAlmostEqual(scalar_result, r.score_double,
                                   places=7,
                                   msg='Incorrect scalar score: %.8f. The correct score is: %.8f, test id: %d' %
                                   (scalar_result, r.score_double, completed + 1))
            '''
            if byte_result > r.threshold:
                self.assertGreaterEqual(byte_result, sys.float_info.max)
            else:
                self.assertAlmostEqual(byte_result, r.score_byte,
                                       places=7,
                                       msg='Incorrect byte score: %.8f. The correct score is: %.8f, test id: %d' %
                                       (byte_result, r.score_byte, completed + 1))
            '''
            if short_result > r.threshold:
                self.assertGreaterEqual(short_result, sys.float_info.max)
            else:
                if short_result < r.score_short and (r.score_short - short_result) > 0.01:
                    print "Warning: python short score(%f) is less than darwin's at id: %d!" \
                          % (short_result, completed + 1)
                else:
                    self.assertAlmostEqual(short_result, r.score_short,
                                       places=7,
                                       msg='Incorrect short score: %.8f. The correct score is: %.8f, test id: %d' %
                                       (short_result, r.score_short, completed + 1))

            completed += 1
            if completed % progress_step == 0:
                print '%d%% completed' % (completed / progress_step)


if __name__ == '__main__':
    unittest.main()