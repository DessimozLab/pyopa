import cython_swps3
import threading
import unittest
import os


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.s1 = 'PDVRTQYSRTKTIKLAQVRKCGAWRVLCLDLIPDLTAKNNHMRTKWTEVQYLAFVVSIVKKRPLSHSLVLITTGKAWNGTWRALPRLSNKLIETAFKEIQAEETVYDTKAFVAGKKPRWVSPFICYGLPFVISRFDFAQYRLKDMLILFSDMLLSRICNFYNGNTGPVPNSKTNEDTDLFFDGLSGMLKLNLKRSDAICHVICYEAPIARVKFGREVKDKFSLPKGGKNPSRRISWNILGILIDRTMFIRPRLVARKEAIHLFDLIGENIDAITQRLRAHKTLMVHESQVVEQPLKVKNLDLRPELVGEEEKNRHGRAKQLDRMANGNMAQIKNGHFKQTYLISVFRPQWLQLQGGCLIAEGFHSEVGGTVDGLKGTPCAQGPVVKGLFAVWRRCDTLAGRYYQKAADIDKLGDILLASLYYIPQGAIITLSEEMAKRIGANVLLVGLINVRYSGIGYEACVGDLAPEVSWLNAGHGNIQMVLHTIDGDGCQTPHGLKIYTDKRLLDLYQGAQLKVTVATTGSVKVSKSMGWLQEGGLDYFALAGRFYRADLREIEHPRAMAVSAHLCAVGLNWVFLADIICDPNEAFKFGKDFEPRTLTYGFANEDENPKNGGATTTSFAVAVYKIKTVATLKVIGKALWKGIQMRTQQGSGPTCQWALRKGKNSILLLAQDSRGGIPKNEFTILGDLPEGQTTTCTHTEIKTRLLYGATVFFMRGDLVGLYADGCSHLYRSSNLMSQACAAAKTILCSLDGERANFSNPTDFAMYNAVFRPRLYTVSFGVFDNNVDVLQAALYYLIMMAMKQYWGVKQGGLEGTLYTWSKVSGKKETSDSRNNPSICVSVCKNPLKDVQLRIAALKRFAEAEEIGKPAVVIRALEPGLTLYILLSSHGSEGKKTHNPILVSAFVVTTVADTSKPKVTYHKDQEMAIYQVLGNNPAGYEVELAFLLPTASSKQQSGRTRKFMDTASGELKEMPIQSSHEITQAADINNLRQLPRTYKKESAKVKVAACKQPPAALNTGIEKVPSHPDGLQLIIEDEWKLLEASSMSQYNEQAKEWPFHKGGIFFKGHEQKCIDASELPRGITRDLRVILINEALVLNTFCGERKLQNEATLILLRAYVWGRHLLANYFRAPNEQDGVLVDIPQGRSTLKSDHLRASIPLFLYTTIETCTSNVTIHKRVQPMIILDIAVAGEGVCDMKNGQVFKRRMARSNDRRLPPGARMKIILFRRNHECYPLQKHQEQWILGAIRTPYGLYNLQEKATLTTRYLIKLQINNRNDLVTTLVSLLMHTRESYIRFTKERRTTESPIDVLAATLYQEFTREVRRAGEQRAGIFFSQDTNYEQAIFETKMAAYPPFGANSWNPTLRYEAWTIIKTPNSKGQEFFLEHMQDVGYGKIASSKYQEKDDDEEVARGRIVPAWY'
        self.s2 = 'PPFQPDKKLAGIELVLCNADLPGRSIYLRKVLQANANKRASASKRCTDDDIIKVDSAPDPQRKLVQAGKVPRVLYNGDVSNIISQILICAYVTGASRNFQHVMLLMDKGWGRGFTLMVNYPCPKVLEEFNPTLLTALVIISVYLNSIECERAGVTIAALNVKLEATDRLALLGRQTANTVMRAPLLLLCQGDSAKNTLNWSLEDLAIVFGRAATRVCKNLALLLNSQVFFQKTTGYKSQLGKNVINFDLYKPLVCDLVDATKYMKFYGTNDDSTDIQGRSSEKAAALAAAAMGVVGWHFLAPTGLVGAGSTFSPVFCIKGNAQLCCKRFDIDEWKALLTLQKSKIANIDYLRYRTGAVIEVGANYDGCSGQPKLQCFYDYLIRYPETVLGTNRQERVMTDEGGEHVRDLILRNVLENPTGFIGSGTHPGNISCTLETTNADLIIGSTDYDGVGSYLIIMGTCFMVTGCVVFTYAVMELVRPLKIHIFACAKVILQEADGSQKTNLRGRGKVSSFGDLPVRFRTLDGIATPSTTHAEMGASFDAAVLVIGRTGTAKFRQFATLDNRNLACNINLSSIRRYFNDNNWLEAGAKNAAEILVNHADKSLTPWVVGLGPLLKPGDIACPMIAVSYLVLVIMDMYLASYSDSFAKHLKNKHRTTTSAHKPSNQQLALDGALTAKRSSQAASIIFEAEEWGFLEWAMIGHLQTKMIYDDAFRLNSPEEELLTQATTHKIKPNYLIALQMLHRDFCIGFFHTLIHASVADSIVYASRLKQNAAIIDRGKTARQDLLGIALKLIVSASTKNAASFNRDFKLPVDVMFRFLDKMLNHGVNTIVHGGQDPKNGNPVGAGLPSWAKNIKVELQVTMFQLFESVDCTSELRLLSTAVDTTLHGEVQVMSAKDLFGRFRYRILSAGESLMENGISPKSFVEALKYFIMYYWTDITEPRCRGSALYPITIQPNLYKRTSATSLHPKGERWLPFEETSRTTISTVLMNNALLGICLYKSYQLLDHDFLGDKKQSNKRVSENSFLGIQTLHDPTGYLQKLDHSRLSKFNRDIRWGQGKSPEQWAVTLVPTLFVKKGTNAWRKKNNAEPIIVTTGTNTAPLEELHKAWMQLAHDGIVVSTLTENEKLEFFSFQDGMPSLVLFSIMAETNQLRYIGNKIYASRKWMADAQKASWVYASLPTNSCNWTAVEVAFEPKGECQMAKKFDLHSMAIVMVRLLAQERSDGADGMNNASSVKWLRKEANEKVCKWWFASPKINAMFQTVKIQSSGKYLARNPKAATKDVKKVEQDLLSRIQTQEHGLLWFYVRLIGEISEVPILSCNKALFLTIKLFNKFIRWNIAPLEITSGVDAWHTIFTSSRFSETDTGIEMTALDLTLPQGNWGTMKKKVALAATGFILFLAYSMGTLSKKFEGNHHWTWVYPFFITITVQLYIFNGHTAWVLFNFVEIPGEAIVSLRTGYLNGGRDKTFVEGLVFNSDVGRTYGGYTSNIK'
        self.s1_norm = cython_swps3.normalize_sequence(self.s1)
        self.s2_norm = cython_swps3.normalize_sequence(self.s2)
        self.envs = cython_swps3.read_all_env_json(
            os.path.dirname(__file__) + '/data/matrices/json/all_matrices.json')
        self.env = self.envs[515]

    def test_norm_seq(self):
        s1 = cython_swps3.normalize_sequence('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        ref = str(bytearray(range(26)))
        self.assertEquals(s1, ref)

        #do not normalize '_' character
        s1 = cython_swps3.normalize_sequence('__________________________', True)
        self.assertEquals(s1, '__________________________')

        self.assertRaises(Exception, cython_swps3.normalize_sequence, 'ABCDEFGHIJLM NOPQRSTUVWXYZ')

    def test_scale(self):
        self.assertEquals(cython_swps3.scale_to_byte(3, 5), 15)
        self.assertEquals(cython_swps3.scale_to_short(3, 5), 15)

        #underflow can occur, an overflow should cause an exception
        self.assertEquals(cython_swps3.scale_to_byte(3, -50), -128)
        self.assertEquals(cython_swps3.scale_to_short(3, -15000), -32768)

        self.assertRaises(Exception, cython_swps3.scale_to_byte, 3, 50)
        self.assertRaises(Exception, cython_swps3.scale_to_short, 3, 15000)

    def test_create_env(self):
        simple_score = 3.78
        simple_scores = [[simple_score]]
        env_a = cython_swps3.create_environment(-2, -1, 20, simple_scores, 'A')
        env_b = cython_swps3.create_environment(-2, -1, 20, simple_scores, 'B')

        self.assertEquals(cython_swps3.align_double('AAA', 'AAA', env_a)[0], 3 * simple_score)
        self.assertEquals(cython_swps3.align_double('BBB', 'BBB', env_b)[0], 3 * simple_score)

        self.assertRaises(Exception, cython_swps3.create_environment, -2, -1, 20, simple_scores, 'AB')
        self.assertRaises(Exception, cython_swps3.create_environment, -2, -1, 20, simple_scores, '')

        bad_matrix = [[0.23], [0.65, -12.32]]
        self.assertRaises(Exception, cython_swps3.create_environment, -2, -1, 20, bad_matrix, 'A')
        self.assertRaises(Exception, cython_swps3.create_environment, -2, -1, 20, bad_matrix, 'AB')

    def test_align(self):
        profile = cython_swps3.AlignmentProfile()
        profile.create_profiles(self.s1, self.env)
        profile_norm = cython_swps3.AlignmentProfile()
        profile_norm.create_profiles(self.s1_norm, self.env, True)

        res_short = cython_swps3.align_short(self.s1, self.s2, self.env)
        res_short_n = cython_swps3.align_short(self.s1_norm, self.s2_norm, self.env, True)
        res_p_short = profile.align_short(self.s2, self.env)
        res_p_short_n = profile.align_short(self.s2_norm, self.env, True)

        res_byte = cython_swps3.align_byte(self.s1, self.s2, self.env)
        res_byte_n = cython_swps3.align_byte(self.s1_norm, self.s2_norm, self.env, True)
        res_p_byte = profile.align_byte(self.s2, self.env)
        res_p_byte_n = profile.align_byte(self.s2_norm, self.env, True)

        res_double = cython_swps3.align_double(self.s1, self.s2, self.env)[0]
        res_double_n = cython_swps3.align_double(self.s1_norm, self.s2_norm, self.env, True)[0]
        res_ref_double = cython_swps3.align_scalar_reference_local(self.s1, self.s2, self.env)
        res_ref_double_n = cython_swps3.align_scalar_reference_local(self.s1_norm, self.s2_norm, self.env, True)

        self.assertAlmostEqual(res_short, res_short_n)
        self.assertAlmostEqual(res_p_short, res_p_short_n)
        self.assertAlmostEqual(res_p_short, res_short)
        self.assertAlmostEqual(res_byte, res_byte_n)
        self.assertAlmostEqual(res_p_byte, res_p_byte_n)
        self.assertAlmostEqual(res_p_byte, res_byte)
        self.assertAlmostEqual(res_double, res_double_n)
        self.assertAlmostEqual(res_ref_double, res_ref_double_n)
        self.assertAlmostEqual(res_ref_double, res_double_n)

        #empty environment should use zero matrix
        self.assertAlmostEqual(cython_swps3.align_double(self.s1, self.s2, cython_swps3.AlignmentEnvironment())[0], 0.0)

    def test_align_strings(self):
        threading.stack_size(100000000)
        t = threading.Thread(None, self._al_strs, 'Align Thread')
        t.start()
        t.join()

    def _al_strs(self):
        alignment_only_max_ranges = cython_swps3.align_double(self.s1, self.s2, self.env, False, False, False, False)
        alignment_full_ranges = cython_swps3.align_double(self.s1, self.s2, self.env, False, False, False, True)

        self.assertRaises(Exception, cython_swps3.align_strings, self.s1, self.s2, self.env,
                          False, False, alignment_only_max_ranges)
        aligned_strings = cython_swps3.align_strings(self.s1, self.s2, self.env, False, False, alignment_full_ranges)
        aligned_strings_norm = cython_swps3.align_strings(self.s1_norm, self.s2_norm, self.env, True)
        aligned_strings_simple = cython_swps3.align_strings(self.s1, self.s2, self.env)

        self.assertEquals(aligned_strings, aligned_strings_norm)
        self.assertEquals(aligned_strings, aligned_strings_simple)

        env = cython_swps3.AlignmentEnvironment()

    def test_estimate_pam(self):
        log_pam1 = cython_swps3.read_env_json(
            os.path.dirname(__file__) + '/data/matrices/json/logPAM1.json')
        dms = cython_swps3.MutipleAlEnv(self.envs, log_pam1)

        self.assertRaises(Exception, dms.estimate_pam, self.s1, self.s2)