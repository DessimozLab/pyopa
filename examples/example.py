import pyopa
import os
import threading


#to do the concrete alignment in a new thread
def nt_align(s1, s2, env, is_global):
    print('Concrete %s alignment:' % ('global' if is_global else 'local'))
    aligned_strings = pyopa.align_strings(s1, s2, env, is_global)
    print('\taligned_s1: %s' % aligned_strings[0])
    print('\taligned_s2: %s' % aligned_strings[1])


def nt_epam(s1, s2, dms, env):
    aligned_strings = pyopa.align_strings(s1, s2, env)
    print('EstimatePam:')
    epam_res = dms.estimate_pam(aligned_strings[0], aligned_strings[1])
    print('\tSim: %f' % epam_res[0])
    print('\tPam Number: %f' % epam_res[1])
    print('\tVariance: %f' % epam_res[2])


s1 = pyopa.Sequence('PDVRTQYSRTKTIKLAQVRKCGAWRVLCLDLIPDLTAKNNHMRTKWTEVQYLAFVVSIVKKRPLSHSLVLITTGKAWNGTWRALPRLSNKLIETAFKEIQAEETVYDTKAFVAGKKPRWVSPFICYGLPFVISRFDFAQYRLKDMLILFSDMLLSRICNFYNGNTGPVPNSKTNEDTDLFFDGLSGMLKLNLKRSDAICHVICYEAPIARVKFGREVKDKFSLPKGGKNPSRRISWNILGILIDRTMFIRPRLVARKEAIHLFDLIGENIDAITQRLRAHKTLMVHESQVVEQPLKVKNLDLRPELVGEEEKNRHGRAKQLDRMANGNMAQIKNGHFKQTYLISVFRPQWLQLQGGCLIAEGFHSEVGGTVDGLKGTPCAQGPVVKGLFAVWRRCDTLAGRYYQKAADIDKLGDILLASLYYIPQGAIITLSEEMAKRIGANVLLVGLINVRYSGIGYEACVGDLAPEVSWLNAGHGNIQMVLHTIDGDGCQTPHGLKIYTDKRLLDLYQGAQLKVTVATTGSVKVSKSMGWLQEGGLDYFALAGRFYRADLREIEHPRAMAVSAHLCAVGLNWVFLADIICDPNEAFKFGKDFEPRTLTYGFANEDENPKNGGATTTSFAVAVYKIKTVATLKVIGKALWKGIQMRTQQGSGPTCQWALRKGKNSILLLAQDSRGGIPKNEFTILGDLPEGQTTTCTHTEIKTRLLYGATVFFMRGDLVGLYADGCSHLYRSSNLMSQACAAAKTILCSLDGERANFSNPTDFAMYNAVFRPRLYTVSFGVFDNNVDVLQAALYYLIMMAMKQYWGVKQGGLEGTLYTWSKVSGKKETSDSRNNPSICVSVCKNPLKDVQLRIAALKRFAEAEEIGKPAVVIRALEPGLTLYILLSSHGSEGKKTHNPILVSAFVVTTVADTSKPKVTYHKDQEMAIYQVLGNNPAGYEVELAFLLPTASSKQQSGRTRKFMDTASGELKEMPIQSSHEITQAADINNLRQLPRTYKKESAKVKVAACKQPPAALNTGIEKVPSHPDGLQLIIEDEWKLLEASSMSQYNEQAKEWPFHKGGIFFKGHEQKCIDASELPRGITRDLRVILINEALVLNTFCGERKLQNEATLILLRAYVWGRHLLANYFRAPNEQDGVLVDIPQGRSTLKSDHLRASIPLFLYTTIETCTSNVTIHKRVQPMIILDIAVAGEGVCDMKNGQVFKRRMARSNDRRLPPGARMKIILFRRNHECYPLQKHQEQWILGAIRTPYGLYNLQEKATLTTRYLIKLQINNRNDLVTTLVSLLMHTRESYIRFTKERRTTESPIDVLAATLYQEFTREVRRAGEQRAGIFFSQDTNYEQAIFETKMAAYPPFGANSWNPTLRYEAWTIIKTPNSKGQEFFLEHMQDVGYGKIASSKYQEKDDDEEVARGRIVPAWY')
s2 = pyopa.Sequence('PPFQPDKKLAGIELVLCNADLPGRSIYLRKVLQANANKRASASKRCTDDDIIKVDSAPDPQRKLVQAGKVPRVLYNGDVSNIISQILICAYVTGASRNFQHVMLLMDKGWGRGFTLMVNYPCPKVLEEFNPTLLTALVIISVYLNSIECERAGVTIAALNVKLEATDRLALLGRQTANTVMRAPLLLLCQGDSAKNTLNWSLEDLAIVFGRAATRVCKNLALLLNSQVFFQKTTGYKSQLGKNVINFDLYKPLVCDLVDATKYMKFYGTNDDSTDIQGRSSEKAAALAAAAMGVVGWHFLAPTGLVGAGSTFSPVFCIKGNAQLCCKRFDIDEWKALLTLQKSKIANIDYLRYRTGAVIEVGANYDGCSGQPKLQCFYDYLIRYPETVLGTNRQERVMTDEGGEHVRDLILRNVLENPTGFIGSGTHPGNISCTLETTNADLIIGSTDYDGVGSYLIIMGTCFMVTGCVVFTYAVMELVRPLKIHIFACAKVILQEADGSQKTNLRGRGKVSSFGDLPVRFRTLDGIATPSTTHAEMGASFDAAVLVIGRTGTAKFRQFATLDNRNLACNINLSSIRRYFNDNNWLEAGAKNAAEILVNHADKSLTPWVVGLGPLLKPGDIACPMIAVSYLVLVIMDMYLASYSDSFAKHLKNKHRTTTSAHKPSNQQLALDGALTAKRSSQAASIIFEAEEWGFLEWAMIGHLQTKMIYDDAFRLNSPEEELLTQATTHKIKPNYLIALQMLHRDFCIGFFHTLIHASVADSIVYASRLKQNAAIIDRGKTARQDLLGIALKLIVSASTKNAASFNRDFKLPVDVMFRFLDKMLNHGVNTIVHGGQDPKNGNPVGAGLPSWAKNIKVELQVTMFQLFESVDCTSELRLLSTAVDTTLHGEVQVMSAKDLFGRFRYRILSAGESLMENGISPKSFVEALKYFIMYYWTDITEPRCRGSALYPITIQPNLYKRTSATSLHPKGERWLPFEETSRTTISTVLMNNALLGICLYKSYQLLDHDFLGDKKQSNKRVSENSFLGIQTLHDPTGYLQKLDHSRLSKFNRDIRWGQGKSPEQWAVTLVPTLFVKKGTNAWRKKNNAEPIIVTTGTNTAPLEELHKAWMQLAHDGIVVSTLTENEKLEFFSFQDGMPSLVLFSIMAETNQLRYIGNKIYASRKWMADAQKASWVYASLPTNSCNWTAVEVAFEPKGECQMAKKFDLHSMAIVMVRLLAQERSDGADGMNNASSVKWLRKEANEKVCKWWFASPKINAMFQTVKIQSSGKYLARNPKAATKDVKKVEQDLLSRIQTQEHGLLWFYVRLIGEISEVPILSCNKALFLTIKLFNKFIRWNIAPLEITSGVDAWHTIFTSSRFSETDTGIEMTALDLTLPQGNWGTMKKKVALAATGFILFLAYSMGTLSKKFEGNHHWTWVYPFFITITVQLYIFNGHTAWVLFNFVEIPGEAIVSLRTGYLNGGRDKTFVEGLVFNSDVGRTYGGYTSNIK')

#loading the matrices and gap costs from JSON
defaults = pyopa.load_default_environments()
envs = defaults['environments']
env = envs[515]

print('Aligning\n%s\nto\n%s\n' % (s1, s2))

#calculating local and global scores for the given sequences
local_double = pyopa.align_double(s1, s2, env)
global_double = pyopa.align_double(s1, s2, env, False, True, True)

#the first element is the score, the other elements of the returned list contain the ranges for the local alignment
print('Local score: %f' % local_double[0])
print('Global score: %f' % global_double[0])

#the align_double function is an efficient vectorized C implementation, however, it is possible to call the
#  reference implementation, and compare the double score given by it to the vectorized version (the scores of course
#  should always be the same)
print('Reference local double score: %f' % pyopa.align_scalar_reference_local(s1, s2, env))

#for the concrete alignment we should increase the stack size
#on linux we can do it by using
#  'resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))'
#or we can start a new thread with the given stack size and do the calculation there
threading.stack_size(100000000)
t = threading.Thread(None, nt_align, 'Aligning Thread', (s1, s2, env, False))
t.start()
t.join()
t = threading.Thread(None, nt_align, 'Aligning Thread', (s1, s2, env, True))
t.start()
t.join()

#we can give an upper bound estimation of the score by using the short and byte alignment
#these estimations are much faster than calculating the actual double score
#for this, we create a profile, which we can use for later estimations
profile = pyopa.AlignmentProfile()
profile.create_profiles(s1, env)
print('Local short estimation with threshold %f: %f' % (env.threshold, profile.align_short(s2, env)))
print('Local byte estimation with threshold %f: %f' % (env.threshold, profile.align_byte(s2, env)))
#since the byte estimation is larger than the threshold, double_max is returned, which indicates that by
#  the byte estimation the score is larger than the threshold
#if we wish to change the threshold, we have to recreate the matrices and the profile used for the estimations,
#  as follows:
env.threshold = 62.0
#if we forget about this part we will get inaccurate results!
env.create_scaled_matrices()
profile.create_profiles(s1, env)
print('Local short estimation with threshold of %f: %f' % (env.threshold, profile.align_short(s2, env)))
print('Local byte estimation with threshold of %f: %f' % (env.threshold, profile.align_byte(s2, env)))

#to hide the profile creation we can also use the simple align_short/byte function, but since
# it always creates a new profile, it is inefficient
print('Simple short alignment %f' % pyopa.align_short(s1, s2, env))
print('Simple byte alignment %f' % pyopa.align_byte(s1, s2, env))

#to use the EstimatePam function we have to create a data structure that stores the matrices

#loading the logPAM1 matrix, which was used to generate the other matrices
log_pam1 = defaults['log_pam1']

dms = pyopa.MutipleAlEnv(envs, log_pam1)
t = threading.Thread(None, nt_epam, 'EstimatePam Thread', (s1, s2, dms, env))
t.start()
t.join()

#to generate environments instead of reading from a file
generated_envs = pyopa.generate_all_env(log_pam1, 1000)

