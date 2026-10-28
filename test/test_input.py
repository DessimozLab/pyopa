import pyopa
import itertools
import numpy
import threading
import unittest
import os
import array
import sys


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.s1 = pyopa.Sequence('PDVRTQYSRTKTIKLAQVRKCGAWRVLCLDLIPDLTAKNNHMRTKWTEVQYLAFVVSIVKKRPLSHSLVLITTGKAWNGTWRALPRLSNKLIETAFKEIQAEETVYDTKAFVAGKKPRWVSPFICYGLPFVISRFDFAQYRLKDMLILFSDMLLSRICNFYNGNTGPVPNSKTNEDTDLFFDGLSGMLKLNLKRSDAICHVICYEAPIARVKFGREVKDKFSLPKGGKNPSRRISWNILGILIDRTMFIRPRLVARKEAIHLFDLIGENIDAITQRLRAHKTLMVHESQVVEQPLKVKNLDLRPELVGEEEKNRHGRAKQLDRMANGNMAQIKNGHFKQTYLISVFRPQWLQLQGGCLIAEGFHSEVGGTVDGLKGTPCAQGPVVKGLFAVWRRCDTLAGRYYQKAADIDKLGDILLASLYYIPQGAIITLSEEMAKRIGANVLLVGLINVRYSGIGYEACVGDLAPEVSWLNAGHGNIQMVLHTIDGDGCQTPHGLKIYTDKRLLDLYQGAQLKVTVATTGSVKVSKSMGWLQEGGLDYFALAGRFYRADLREIEHPRAMAVSAHLCAVGLNWVFLADIICDPNEAFKFGKDFEPRTLTYGFANEDENPKNGGATTTSFAVAVYKIKTVATLKVIGKALWKGIQMRTQQGSGPTCQWALRKGKNSILLLAQDSRGGIPKNEFTILGDLPEGQTTTCTHTEIKTRLLYGATVFFMRGDLVGLYADGCSHLYRSSNLMSQACAAAKTILCSLDGERANFSNPTDFAMYNAVFRPRLYTVSFGVFDNNVDVLQAALYYLIMMAMKQYWGVKQGGLEGTLYTWSKVSGKKETSDSRNNPSICVSVCKNPLKDVQLRIAALKRFAEAEEIGKPAVVIRALEPGLTLYILLSSHGSEGKKTHNPILVSAFVVTTVADTSKPKVTYHKDQEMAIYQVLGNNPAGYEVELAFLLPTASSKQQSGRTRKFMDTASGELKEMPIQSSHEITQAADINNLRQLPRTYKKESAKVKVAACKQPPAALNTGIEKVPSHPDGLQLIIEDEWKLLEASSMSQYNEQAKEWPFHKGGIFFKGHEQKCIDASELPRGITRDLRVILINEALVLNTFCGERKLQNEATLILLRAYVWGRHLLANYFRAPNEQDGVLVDIPQGRSTLKSDHLRASIPLFLYTTIETCTSNVTIHKRVQPMIILDIAVAGEGVCDMKNGQVFKRRMARSNDRRLPPGARMKIILFRRNHECYPLQKHQEQWILGAIRTPYGLYNLQEKATLTTRYLIKLQINNRNDLVTTLVSLLMHTRESYIRFTKERRTTESPIDVLAATLYQEFTREVRRAGEQRAGIFFSQDTNYEQAIFETKMAAYPPFGANSWNPTLRYEAWTIIKTPNSKGQEFFLEHMQDVGYGKIASSKYQEKDDDEEVARGRIVPAWY')
        self.s2 = pyopa.Sequence('PPFQPDKKLAGIELVLCNADLPGRSIYLRKVLQANANKRASASKRCTDDDIIKVDSAPDPQRKLVQAGKVPRVLYNGDVSNIISQILICAYVTGASRNFQHVMLLMDKGWGRGFTLMVNYPCPKVLEEFNPTLLTALVIISVYLNSIECERAGVTIAALNVKLEATDRLALLGRQTANTVMRAPLLLLCQGDSAKNTLNWSLEDLAIVFGRAATRVCKNLALLLNSQVFFQKTTGYKSQLGKNVINFDLYKPLVCDLVDATKYMKFYGTNDDSTDIQGRSSEKAAALAAAAMGVVGWHFLAPTGLVGAGSTFSPVFCIKGNAQLCCKRFDIDEWKALLTLQKSKIANIDYLRYRTGAVIEVGANYDGCSGQPKLQCFYDYLIRYPETVLGTNRQERVMTDEGGEHVRDLILRNVLENPTGFIGSGTHPGNISCTLETTNADLIIGSTDYDGVGSYLIIMGTCFMVTGCVVFTYAVMELVRPLKIHIFACAKVILQEADGSQKTNLRGRGKVSSFGDLPVRFRTLDGIATPSTTHAEMGASFDAAVLVIGRTGTAKFRQFATLDNRNLACNINLSSIRRYFNDNNWLEAGAKNAAEILVNHADKSLTPWVVGLGPLLKPGDIACPMIAVSYLVLVIMDMYLASYSDSFAKHLKNKHRTTTSAHKPSNQQLALDGALTAKRSSQAASIIFEAEEWGFLEWAMIGHLQTKMIYDDAFRLNSPEEELLTQATTHKIKPNYLIALQMLHRDFCIGFFHTLIHASVADSIVYASRLKQNAAIIDRGKTARQDLLGIALKLIVSASTKNAASFNRDFKLPVDVMFRFLDKMLNHGVNTIVHGGQDPKNGNPVGAGLPSWAKNIKVELQVTMFQLFESVDCTSELRLLSTAVDTTLHGEVQVMSAKDLFGRFRYRILSAGESLMENGISPKSFVEALKYFIMYYWTDITEPRCRGSALYPITIQPNLYKRTSATSLHPKGERWLPFEETSRTTISTVLMNNALLGICLYKSYQLLDHDFLGDKKQSNKRVSENSFLGIQTLHDPTGYLQKLDHSRLSKFNRDIRWGQGKSPEQWAVTLVPTLFVKKGTNAWRKKNNAEPIIVTTGTNTAPLEELHKAWMQLAHDGIVVSTLTENEKLEFFSFQDGMPSLVLFSIMAETNQLRYIGNKIYASRKWMADAQKASWVYASLPTNSCNWTAVEVAFEPKGECQMAKKFDLHSMAIVMVRLLAQERSDGADGMNNASSVKWLRKEANEKVCKWWFASPKINAMFQTVKIQSSGKYLARNPKAATKDVKKVEQDLLSRIQTQEHGLLWFYVRLIGEISEVPILSCNKALFLTIKLFNKFIRWNIAPLEITSGVDAWHTIFTSSRFSETDTGIEMTALDLTLPQGNWGTMKKKVALAATGFILFLAYSMGTLSKKFEGNHHWTWVYPFFITITVQLYIFNGHTAWVLFNFVEIPGEAIVSLRTGYLNGGRDKTFVEGLVFNSDVGRTYGGYTSNIK')

        defaults = pyopa.load_default_environments()
        self.envs = defaults['environments']
        self.log_pam1 = defaults['log_pam1']
        self.env = self.envs[515]

    def test_norm_seq(self):
        s1 = pyopa.normalize_sequence('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        if (sys.version_info < (3,)):
            ref = bytearray(range(26))
            self.assertEqual(s1, ref)
        else:
            ref = bytes(range(26)).decode('utf-8')
            self.assertEqual(s1, ref)

        # do not normalize '_' character
        s1 = pyopa.normalize_sequence('__________________________')
        self.assertEqual(s1, '__________________________')

        self.assertRaises(ValueError, pyopa.normalize_sequence, 'ABCDEFGHIJLM NOPQRSTUVWXYZ')

    def test_scale(self):
        self.assertEqual(pyopa.scale_to_byte(3, 5), 15)
        self.assertEqual(pyopa.scale_to_short(3, 5), 15)

        # underflow can occur, an overflow should cause an exception
        self.assertEqual(pyopa.scale_to_byte(3, -50), -128)
        self.assertEqual(pyopa.scale_to_short(3, -15000), -32768)

        self.assertRaises(OverflowError, pyopa.scale_to_byte, 3, 50)
        self.assertRaises(OverflowError, pyopa.scale_to_short, 3, 15000)

    def test_scaled_matrices(self):
        # the vectorized scaling must match the scalar functions entry by entry
        for env in self.envs[::25]:
            for threshold in (40.0, 85.0, 300.0):
                env.threshold = threshold
                env.create_scaled_matrices()
                byte_factor, short_factor = env.byte_factor(), env.short_factor()
                for i in range(26):
                    for j in range(26):
                        val = env.float64_matrix[i][j]
                        self.assertEqual(env.int8_matrix[i][j], pyopa.scale_to_byte(val, byte_factor))
                        self.assertEqual(env.int16_matrix[i][j], pyopa.scale_to_short(val, short_factor))
                self.assertEqual(env.int8_matrix.dtype, numpy.int8)
                self.assertEqual(env.int16_matrix.dtype, numpy.int16)

        # a threshold this low makes the factor so large that the matrix overflows
        env = self.envs[515]
        env.threshold = 1.0
        self.assertRaises(OverflowError, env.create_scaled_matrices)
        env.threshold = 85.0

    def test_ranges_match_full_reverse_scan(self):
        # the reverse pass for the ranges stops early once it reaches the forward score;
        # it must find the same start as scanning the whole reversed prefix
        seqs = [self.s1, self.s2, pyopa.Sequence(self.s1.convert_readable()[300:900]),
                pyopa.Sequence(self.s2.convert_readable()[100:400])]
        for env in self.envs[::200]:
            for a, b in itertools.permutations(seqs, 2):
                score, max1, max2, min1, min2 = pyopa.align_double(a, b, env, False, False, True)

                profile = pyopa.AlignmentProfile()
                profile.create_profile_double(pyopa.Sequence(a.s_norm[max1::-1], True), env.float64_matrix)
                rev = profile.align_double(pyopa.Sequence(b.s_norm[max2::-1], True), env)
                self.assertAlmostEqual(rev[0], score, delta=abs(score) * 1e-10)
                self.assertEqual((min1, min2), (max1 - rev[1], max2 - rev[2]))

    def test_ranges_near_tie(self):
        # with A/A almost 0, the alignment without the leading A(s) scores only marginally
        # less than the full one; the reverse pass must not stop at that shorter alignment
        env = pyopa.create_environment(-2, -1, 20, [[1e-11, -1], [-1, 1]], 'AB')
        for s in ['AB', 'ABBBBB', 'AAB', 'A' + 'B' * 40]:
            seq = pyopa.Sequence(s)
            score, max1, max2, min1, min2 = pyopa.align_double(seq, seq, env, False, False, True)
            self.assertEqual((min1, min2, max1, max2), (0, 0, len(s) - 1, len(s) - 1), msg=s)

    def test_create_env(self):
        simple_score = 3.78
        simple_scores = [[simple_score]]
        env_a = pyopa.create_environment(-2, -1, 20, simple_scores, 'A')
        env_b = pyopa.create_environment(-2, -1, 20, simple_scores, 'B')

        s_short1 = pyopa.Sequence('AAA')
        s_short2 = pyopa.Sequence('BBB')

        s_short3 = pyopa.Sequence('A')
        s_short4 = pyopa.Sequence('B')
        s_short5 = pyopa.Sequence('')

        self.assertEqual(pyopa.align_double(s_short1, s_short1, env_a)[0], 3 * simple_score)
        self.assertEqual(pyopa.align_double(s_short2, s_short2, env_b)[0], 3 * simple_score)
        self.assertEqual(pyopa.align_double(s_short3, s_short3, env_a)[0], 1 * simple_score)
        self.assertEqual(pyopa.align_double(s_short4, s_short4, env_b)[0], 1 * simple_score)
        self.assertEqual(pyopa.align_double(s_short3, s_short5, env_a)[0], 0.0)
        self.assertEqual(pyopa.align_double(s_short5, s_short5, env_a)[0], 0.0)
        self.assertEqual(pyopa.align_double(s_short3, s_short1, env_a)[0], 1 * simple_score)
        self.assertEqual(pyopa.align_double(s_short1, s_short3, env_a)[0], 1 * simple_score)

        self.assertRaises(ValueError, pyopa.create_environment, -2, -1, 20, simple_scores, 'AB')
        self.assertRaises(ValueError, pyopa.create_environment, -2, -1, 20, simple_scores, '')

        bad_matrix = [[0.23], [0.65, -12.32]]
        self.assertRaises(ValueError, pyopa.create_environment, -2, -1, 20, bad_matrix, 'A')
        self.assertRaises(ValueError, pyopa.create_environment, -2, -1, 20, bad_matrix, 'AB')

    def test_align(self):
        profile = pyopa.AlignmentProfile()
        profile.create_profiles(self.s1, self.env)

        res_short = pyopa.align_short(self.s1, self.s2, self.env)
        res_p_short = profile.align_short(self.s2, self.env)

        res_byte = pyopa.align_byte(self.s1, self.s2, self.env)
        res_p_byte = profile.align_byte(self.s2, self.env)

        res_double = pyopa.align_double(self.s1, self.s2, self.env)[0]
        res_ref_double = pyopa.align_scalar_reference_local(self.s1, self.s2, self.env)

        self.assertAlmostEqual(res_p_short, res_short)
        self.assertAlmostEqual(res_p_byte, res_byte)
        self.assertAlmostEqual(res_double, res_ref_double)

        #empty environment should use zero matrix
        self.assertAlmostEqual(pyopa.align_double(self.s1, self.s2, pyopa.AlignmentEnvironment())[0], 0.0)
        pyopa.align_double(self.s1, self.s2, self.env, True, True, True)

    #def test_align_strings(self):
    #    threading.stack_size(4096*256*64)
    #    t = threading.Thread(None, self._al_strs, 'Align Thread')
    #    t.start()
    #    t.join()

    #def _al_strs(self):
    def test_align_strings(self):
        alignment_only_max_ranges = pyopa.align_double(self.s1, self.s2, self.env, False, False, False)
        alignment_full_ranges = pyopa.align_double(self.s1, self.s2, self.env, False, False, True)

        self.assertRaises(ValueError, pyopa.align_strings, self.s1, self.s2, self.env,
                          False, alignment_only_max_ranges)
        aligned_strings = pyopa.align_strings(self.s1, self.s2, self.env, False, alignment_full_ranges)
        aligned_strings_norm = pyopa.align_strings(self.s1, self.s2, self.env)

        self.assertEqual(aligned_strings, aligned_strings_norm)

        # check __ne__
        self.assertNotEqual(aligned_strings[0], aligned_strings[1])

    def test_estimate_pam(self):
        dms = pyopa.MutipleAlEnv(self.envs, self.log_pam1)

        self.assertRaises(ValueError, dms.estimate_pam, self.s1, self.s2)

    def test_sequence(self):
        s_string = 'TE_ST'
        s_normalized = pyopa.normalize_sequence(s_string)

        # checking non-normalized constructor
        s = pyopa.Sequence(s_string)
        self.assertEqual(s.convert_readable(), s_string)

        # normalized constructor
        s = pyopa.Sequence(s_normalized, True)
        self.assertEqual(s.convert_readable(), s_string)

        # wrong type exception
        if (sys.version_info < (3,)): 
            s_bytes = array.array('B', s_string)
            self.assertRaises(ValueError, pyopa.Sequence, s_bytes, True)
        else:
            s_bytes = s_string.encode('utf-8')
            self.assertRaises(ValueError, pyopa.Sequence, s_bytes, True)
            
        # normalized and non-normalized byte list constructor
        self.assertEqual('ACA_', pyopa.Sequence([0, 2, 0, ord('_')], True))
        self.assertEqual('ACA_', pyopa.Sequence([65, 67, 65, ord('_')], False))


if __name__ == '__main__':
    unittest.main()
