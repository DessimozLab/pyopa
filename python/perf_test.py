import cython_swps3
import time
import profile

def test_func():
    env = cython_swps3.readEnvJson('516.dat')

    with open('test/seq.test') as f:
        sequences = f.read().splitlines()

    sequences = map(lambda s: cython_swps3.normalizeString(s), sequences)

    result = []


    start_time = time.time()

    for i in range(len(sequences)):
        profile = cython_swps3.AlignmentProfile()
        profile.createProfilesNormalized(sequences[i], env)

        for j in range(i + 1, len(sequences)):
            result.append(profile.alignShortNormalized(sequences[j], env))


    elapsed = time.time() - start_time

    print "A total of %d alignments have been done in %.3fs" % (len(result), elapsed)


#profile.run('test_func(); print')
test_func()