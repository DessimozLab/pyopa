import unittest
import pyopa
import numpy
import os


class EnvironmentCreatorTest(unittest.TestCase):
    
    def setUp(self):
       self.log_pam1 = pyopa.read_env_json(os.path.join(pyopa.matrix_dir(), 'logPAM1.json'))
       self.env_generated = pyopa.generate_all_env(self.log_pam1, 1266)
       self.env_loaded = pyopa.load_default_environments()['environments']

    def test_equal_envs(self):
        self.assertEqual(len(self.env_generated), len(self.env_loaded))
        for i in range(len(self.env_loaded)):
            self.assertTrue(numpy.allclose(self.env_generated[i].float64_matrix,
                                              self.env_loaded[i].float64_matrix, 
                                              rtol=1e-10, atol=1e-9),
                            'matrix[{}] differs between loaded and generated:\n {}'
                            .format(i, numpy.max(abs(self.env_generated[i].float64_matrix - 
                                self.env_loaded[i].float64_matrix))))
            for memb in ['threshold', 'gap_open', 'gap_ext', 'pam']:
                attr_gen = getattr(self.env_generated[i], memb)
                attr_load = getattr(self.env_loaded[i], memb)
                self.assertAlmostEqual(attr_gen, attr_load, delta=1e-10*abs(attr_load),
                                       msg='attr {} differ for env[{}]: {:.12g} != {:.12g}'
                                           .format(memb, i, attr_gen, attr_load))


if __name__=='__main__':
    unittest.main()
