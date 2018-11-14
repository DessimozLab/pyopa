import unittest
import os
import csv
import pyopa
import sys
import json
import threading
import numpy as np

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
        self.als1 = ''
        self.als2 = ''
        self.ep_sim = 0
        self.ep_pamn = 0
        self.ep_var = 0

"""
def write_env_file(env, name):
    matrix = env.float64_matrix.reshape(26*26)
    with open(os.path.dirname(__file__) + '/data/matrices/C_compatible/' + name + '.dat', 'w+') as f:
        f.write("%.20f\n%.20f\n%.20f\n" % (env.gap_open, env.gap_ext, env.pam))

        for i in range(26*26):
            f.write("%.20f " % matrix[i])

def write_all_env_files(all_envs):

    for i in range(len(all_envs)):
        write_env_file(all_envs[i], str(i + 1))
"""


class AlignTest(unittest.TestCase):

    def setUp(self):
        self.precision = 10
        #resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))

        with open(os.path.join(os.path.dirname(__file__), 'data', 'testseqs.txt')) as f:
            self.sequences = f.readlines()

        self.sequences = list(map(lambda s: pyopa.Sequence(s.strip()), self.sequences))
        self.darwin_results = []

        defaults = pyopa.load_default_environments()
        self.alignment_environments = defaults['environments']
        self.log_pam1 = defaults['log_pam1']
        self.dms = pyopa.MutipleAlEnv(self.alignment_environments, self.log_pam1)

        """
        write_all_env_files(self.alignment_environments)
        with open(os.path.dirname(__file__) + '/data/matrices/json/logPAM1.json') as lp:
            json_data = json.load(lp)
            #json_data["Scores"] = map(lambda l: map(lambda s: s/(2048*2048*2048), l), json_data["Scores"])
            logPAM1 = pyopa.read_env_json(json_data, self.alignment_environments[0].columns)
            write_env_file(logPAM1, "logPAM1")
        """

        with open(os.path.join(os.path.dirname(__file__), 'data', 'reference_test_results.dat')) as f:
            #skip header
            next(f)
            reader = csv.reader(f, delimiter='\t')

            for s1, s2, matrix_nr, pam, threshold, score_d, score_f, score_s,\
                score_b, als1, als2, ep_sim, ep_pamn, ep_var, in reader:
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
                curr.als1 = als1
                curr.als2 = als2
                curr.ep_sim = float(ep_sim)
                curr.ep_pamn = float(ep_pamn)
                curr.ep_var = float(ep_var)

                self.darwin_results.append(curr)

                '''
                if curr.s1_id not in self.alignment_profiles:
                    p = pyopa.AlignmentProfile()
                    p.create_profiles(self.sequences[curr.s1_id - 1], self.alignment_environments[curr.s1_id - 1])
                    self.alignment_profiles[curr.s1_id] = p
                '''

    #def test_align(self):
    #    threading.stack_size(67104768)
    #    t = threading.Thread(None, self._align_t, 'Aligning Thread')
    #    t.start()
    #    t.join()

    #def _align_t(self):
    @unittest.skip("skip long runningtest")
    def test_align(self):
        print('Running alignment tests...')
        completed = 0
        max_alignments = len(self.darwin_results)
        progress_step = int(np.ceil(max_alignments / 100))

        for r in self.darwin_results:
            s1 = self.sequences[r.s1_id - 1]
            s2 = self.sequences[r.s2_id - 1]
            env = self.alignment_environments[r.matrix_nr - 1]
            env.threshold = r.threshold
            env.create_scaled_matrices()

            #profile = self.alignment_profiles[r.s1_id]
            profile = pyopa.AlignmentProfile()
            profile.create_profiles(s1, env)

            scalar_result_reference = pyopa.align_scalar_reference_local(s1, s2, env)
            double_alignment = pyopa.align_double(s1, s2, env, False, False, True)
            double_result = double_alignment[0]
            byte_result = profile.align_byte(s2, env)
            short_result = profile.align_short(s2, env)

            if r.als1 != '':
                aligned_strings = pyopa.align_strings(s1, s2, env, False, double_alignment)
                ep_result = self.dms.estimate_pam(aligned_strings[0], aligned_strings[1])
                self.assertEqual(aligned_strings[0], r.als1)
                self.assertEqual(aligned_strings[1], r.als2)
                self.assertAlmostEqual(ep_result[0], r.ep_sim,
                                       delta=r.ep_sim *10**(1-self.precision),
                                       msg='Incorrect EstimatePam similarity score: %.10f != %.10f.'
                                           'Test id: %d' %
                                           (ep_result[0], r.ep_sim, completed + 1))
                self.assertAlmostEqual(ep_result[1], r.ep_pamn, delta=r.ep_pamn*10**(1-self.precision))
                self.assertAlmostEqual(ep_result[2], r.ep_var, delta=r.ep_var*10**(1-self.precision))

            self.assertAlmostEqual(scalar_result_reference, r.score_double, places=self.precision,
                                   msg='Incorrect reference double score: %.8f. The correct score is: %.8f, test id: %d'
                                       % (scalar_result_reference, r.score_double, completed + 1))
            self.assertGreaterEqual(short_result, r.score_double,
                                    msg="Short score (%f) must be greater or equal"
                                        " than double score (%f), test id: %d'."
                                        % (short_result, r.score_double, completed + 1)
            )
            self.assertGreaterEqual(byte_result, r.score_double,
                                    msg="Byte score must be greater or equal than double score.")
            self.assertAlmostEqual(double_result, r.score_double,
                                   places=self.precision,
                                   msg='Incorrect double score: %.8f. The correct score is: %.8f, test id: %d' %
                                   (double_result, r.score_double, completed + 1))
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
                    #print("Warning: python short score(%f) is less than darwin's, but still bigger" \
                    #      " than the double score(%f) at id: %d!" \
                    #      % (short_result, r.score_double, completed + 1))
                    pass
                else:
                    self.assertAlmostEqual(short_result, r.score_short,
                                           places=self.precision,
                                           msg='Incorrect short score: %.8f. The correct score is: %.8f, test id: %d' %
                                               (short_result, r.score_short, completed + 1))

            completed += 1
            if (completed % progress_step == 0):
                print('%d%% completed' % (completed / progress_step))

    def test_generated_envs(self):
        print('Testing generated matrices')
        generated_envs = pyopa.generate_all_env(self.log_pam1, 1266)
        for i in range(1266):
            curr_ref = self.alignment_environments[i]
            curr_gen = generated_envs[i]
            self.assertAlmostEqual(curr_ref.pam, curr_gen.pam)
            self.assertAlmostEqual(curr_ref.threshold, curr_gen.threshold)
            self.assertAlmostEqual(curr_ref.gap_open, curr_gen.gap_open)
            self.assertAlmostEqual(curr_ref.gap_ext, curr_gen.gap_ext)

            for j in range(26):
                for k in range(26):
                    self.assertAlmostEqual(curr_ref.float64_matrix[j][k], curr_gen.float64_matrix[j][k])


if __name__ == '__main__':
    unittest.main()
