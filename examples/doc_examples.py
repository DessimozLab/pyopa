import pyopa
import os
import threading


#---------------------------------------------------------------------------------------------------
data = {'gap_open': -20.56,
        'gap_ext': -3.37,
        'pam_distance': 150.87,
        'scores': [[10.0]],
        'column_order': 'A',
        'threshold': 50.0}

env = pyopa.create_environment(**data)

s1 = pyopa.Sequence('AAA')
s2 = pyopa.Sequence('TTT')

#prints [30.0, 2, 2, 0, 0], the first element is the score
print(pyopa.align_double(s1, s1, env))

#prints [0.0, -1, -1, 0, 0], the score is 0
# since the score for 'A -> T' is undefined
print(pyopa.align_double(s2, s1, env))

#---------------------------------------------------------------------------------------------------
#loading the default environments from the data directory
# created at installation time
defaults = pyopa.load_default_environments()
env_list = defaults['environments']
log_pam1_env = defaults['log_pam1']

#the default directory (created at installation time)
matrix_dir = pyopa.matrix_dir()

#or alternatively, you can specify an exact location
env_list = pyopa.read_all_env_json(
    os.path.join(matrix_dir, 'all_matrices.json'))
log_pam1_env = pyopa.read_env_json(
    os.path.join(matrix_dir, 'logPAM1.json'))
#---------------------------------------------------------------------------------------------------
#generates a signle AlignmentEnvironment
# with a pam distance of 250
generated_env = pyopa.generate_env(log_pam1_env, 250)

#generates 1000 environments for different pam distances
gen_env_list = pyopa.generate_all_env(log_pam1_env, 1000)
#---------------------------------------------------------------------------------------------------
s1 = pyopa.Sequence('AATCGGA')
s2 = pyopa.Sequence('AAAA')
s3 = pyopa.Sequence('CATACCTGGTGTGATGCC')

#not optimal, multiple hidden profile generations in the background
print(pyopa.align_short(s1, s2, generated_env))
print(pyopa.align_short(s1, s3, generated_env))
print(pyopa.align_byte(s1, s2, generated_env))
print(pyopa.align_byte(s1, s3, generated_env))

#one profile generation
profile = pyopa.AlignmentProfile()
profile.create_profiles(s1, generated_env)

#the following code produces the exact same result
#but is more efficient since it's
#using the same profile for multiple alignments
print(profile.align_short(s2, generated_env))
print(profile.align_short(s3, generated_env))
print(profile.align_byte(s2, generated_env))
print(profile.align_byte(s3, generated_env))
#---------------------------------------------------------------------------------------------------
generated_env.threshold = 50.0
#WRONG: the short and byte matrices are not recomputed!!!
# the matrices and gap costs for the old threshold
# will be used
profile.create_profiles(s1, generated_env)

#recompute the byte/short matrices and gap costs
generated_env.create_scaled_matrices()

#and then create the profile
profile.create_profiles(s1, generated_env)
#now we can do alignments with the new profile:
print(profile.align_short(s3, generated_env))
#---------------------------------------------------------------------------------------------------
#always a local alignment
print(pyopa.align_scalar_reference_local(s1, s2, generated_env))
#---------------------------------------------------------------------------------------------------
s1_norm = pyopa.normalize_sequence('AATCGGA')

#if the sequence comes from a normalized source
s1 = pyopa.Sequence(s1_norm, True)
s2 = pyopa.Sequence('AAAA')

#construct from a byte array, prints ACCA
print(pyopa.Sequence([0, 2, 2, 0], True))

print(pyopa.align_double(s1, s2, generated_env))
#---------------------------------------------------------------------------------------------------
s1 = pyopa.Sequence('AATCGGA')
s3 = pyopa.Sequence('CATACCTGGTGTGATGCC')
#does not stop at threshold,
#  it is NOT a global alignment, and computes the ranges
#returns [19.946, 6, 8, 0, 2], the first element is the score
# the 4 other elements are [max1, max2, min1, min2] the ranges
print(pyopa.align_double(s1, s3, generated_env, False, False, True))

#returns [score, max1, max2] because the last flag (calculate ranges) is false
print(pyopa.align_double(s1, s3, generated_env, False, False, False))

generated_env.threshold = 10.0
#no generated_env.create_scaled_matrices() is needed
# because we do not operate on short/byte matrices

#results in [11.499268729503227, 3, 0, 3, 0], not the best possible
# local alignment, but still over the threshold of 10.0
print(pyopa.align_double(s1, s3, generated_env, True, False, True))

#global alignment, stop at threshold is ignored
print(pyopa.align_double(s1, s3, generated_env, True, True, True))
#---------------------------------------------------------------------------------------------------
#to do the concrete alignment in a new thread
#or alternatively you can increase your stack size on UNIX-based systems:
#'resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))'
def nt_align(s1, s2, env, is_global, aligned_strs):
    print('Concrete %s alignment:' % ('global' if is_global else 'local'))
    tmp_aligned_strings = pyopa.align_strings(s1, s2, env, is_global)
    print('\taligned_s1: %s' % tmp_aligned_strings[0])
    print('\taligned_s2: %s' % tmp_aligned_strings[1])
    aligned_strs.extend(tmp_aligned_strings)

s1 = pyopa.Sequence('PISRIDNNKITTTLGNTGIISVTIGVIIFKDLHAKVHGF')
s2 = pyopa.Sequence('PIERIENNKILANTGVISVTIGVIIYQDLHADTVMTSDY')
threading.stack_size(100000000)

# aligned_s1: PISRIDNNKITTTLGNTGIISVTIGVIIFKDLHAKV
# aligned_s2: PIERIENNKI___LANTGVISVTIGVIIYQDLHADT
aligned_strings = []
t = threading.Thread(None, nt_align,
                     'Aligning Thread', (s1, s2, generated_env, False, aligned_strings))
t.start()
t.join()
print(aligned_strings[0])
print(aligned_strings[1])
#---------------------------------------------------------------------------------------------------
dms = pyopa.MutipleAlEnv(gen_env_list, log_pam1_env)

#returns an array: [similarity, pam_distance, variance]
print(dms.estimate_pam(aligned_strings[0], aligned_strings[1]))


