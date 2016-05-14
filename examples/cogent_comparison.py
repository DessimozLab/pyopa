import pyopa
import os
import time


def measure_performance(align_func, sequences, env):
    start_time = time.time()
    align_func(sequences, env)
    elapsed_time = time.time() - start_time

    return elapsed_time


def all_against_all_double(sequences, env):
    seq_num = len(sequences)
    for i in range(seq_num):
        profile = pyopa.AlignmentProfile()
        profile.create_profile_double(sequences[i], env.float64_matrix)

        for j in range(i + 1, seq_num):
            profile.align_double(sequences[j], env)


def all_against_all_double_old(sequences, env):
    seq_num = len(sequences)
    for i in range(seq_num):
        s1 = sequences[i]
        for j in range(i + 1, seq_num):
            pyopa.align_double(s1, sequences[j], env, False, False, False)


def all_against_all_byte(sequences, env):
    seq_num = len(sequences)
    for i in range(seq_num):
        profile = pyopa.AlignmentProfile()
        profile.create_profile_byte(sequences[i], env.int8_matrix)

        for j in range(i + 1, seq_num):
            profile.align_byte(sequences[j], env)


def all_against_all_short(sequences, env):
    seq_num = len(sequences)
    for i in range(seq_num):
        profile = pyopa.AlignmentProfile()
        profile.create_profile_short(sequences[i], env.int16_matrix)

        for j in range(i + 1, seq_num):
            profile.align_short(sequences[j], env)


def read_sequences(seq_file):
    with open(seq_file) as f:
        sequences = f.readlines()

    return [pyopa.Sequence(s.strip()) for s in sequences]


seqs = read_sequences(os.path.join(os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../test/data/cmp_seqs.txt')))
env = pyopa.load_default_environments()['environments'][515]


print('All against all (double precision): %fs' % measure_performance(all_against_all_double, seqs, env))
print('All against all (double precision old): %fs' % measure_performance(all_against_all_double_old, seqs, env))
print('All against all (byte estimation): %fs' % measure_performance(all_against_all_byte, seqs, env))
print('All against all (short estimation): %fs' % measure_performance(all_against_all_short, seqs, env))
