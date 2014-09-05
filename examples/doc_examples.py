import cython_swps3
import os


#---------------------------------------------------------------------------------------------------
data = {'gap_open': -20.56,
        'gap_ext': -3.37,
        'pam_distance': 150.87,
        'scores': [[10.0]],
        'column_order': 'A',
        'threshold': 50.0}

env = cython_swps3.create_environment(**data)

#prints [30.0, 2, 2, 0, 0], the first element is the score
print cython_swps3.align_double('AAA', 'AAA', env)

#prints [0.0, -1, -1, 0, 0], the score is 0
# since the score for 'A -> T' is undefined
print cython_swps3.align_double('TTT', 'AAA', env)

#---------------------------------------------------------------------------------------------------
data_dir = os.path.dirname(__file__) +\
           '/../test/data/matrices/json/'
env_list = cython_swps3.read_all_env_json(
    os.path.join(data_dir, 'all_matrices.json'))
log_pam1_env = cython_swps3.read_env_json(
    os.path.join(data_dir, 'logPAM1.json'))
#---------------------------------------------------------------------------------------------------
#generates a signle AlignmentEnvironment
# with a pam distance of 250
generated_env = cython_swps3.generate_env(log_pam1_env, 250)

#generates 1000 environments for different pam distances
gen_env_list = cython_swps3.generate_all_env(log_pam1_env, 1000)
#---------------------------------------------------------------------------------------------------
s1 = 'AATCGGA'
s2 = 'AAAA'
s3 = 'CATACCTGGTGTGATGCC'

#not optimal, two hidden profile generation in the background
print cython_swps3.align_short(s1, s2, generated_env)
print cython_swps3.align_short(s1, s3, generated_env)
print cython_swps3.align_byte(s1, s2, generated_env)
print cython_swps3.align_byte(s1, s3, generated_env)

#one profile generation
profile = cython_swps3.AlignmentProfile()
profile.create_profiles(s1, generated_env)

#the following code produces the exact same result
#but is more efficient since it's
#using the same profile for multiple alignments
print profile.align_short(s2, generated_env)
print profile.align_short(s3, generated_env)
print profile.align_byte(s2, generated_env)
print profile.align_byte(s3, generated_env)
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
print profile.align_short(s3, generated_env)
#---------------------------------------------------------------------------------------------------
#always a local alignment
print cython_swps3.align_scalar_reference_local(s1, s2, generated_env)
#---------------------------------------------------------------------------------------------------
s1_norm = cython_swps3.normalize_sequence(s1)
s2_norm = cython_swps3.normalize_sequence(s2)
s3_norm = cython_swps3.normalize_sequence(s3)

#normalizing sequences in the background
print cython_swps3.align_double(s1, s2, generated_env)

#passing normalized sequences and setting the is_normalized
# flag, a more efficient solution
print cython_swps3.align_double(s1_norm, s2_norm, generated_env, True)
#---------------------------------------------------------------------------------------------------
s1 = 'AATCGGA'
s3 = 'CATACCTGGTGTGATGCC'
#operating on normalized sequences, does not stop at threshold,
#  it is NOT a global alignment, and computes the ranges
#returns [19.946195452221108, 6, 8, 0, 2], the first element is the score
# the 4 other elements are [max1, max2, min1, min2] the ranges
print cython_swps3.align_double(s1_norm, s3_norm, generated_env, True, False, False, True)

#returns [score, max1, max2]
print cython_swps3.align_double(s1_norm, s3_norm, generated_env, True, False, False, False)

generated_env.threshold = 10.0
#no generated_env.create_scaled_matrices() is needed
# because we do not operate on short/byte matrices

#results in [11.499268729503227, 3, 0, 3, 0], not the best possible
# local alignment, but still over the threshold of 10.0
print cython_swps3.align_double(s1_norm, s3_norm, generated_env, True, True, False, True)

#global alignment, stop at threshold is ignored
print cython_swps3.align_double(s1, s3, generated_env, False, True, True, True)
#---------------------------------------------------------------------------------------------------