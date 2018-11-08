import json
import logging
import math
import os
import pkgutil
import re
import sys
from math import log10
from subprocess import call

cimport cpyopa
cimport numpy as np
import numpy as np


def matrix_dir():
    default_matrix_dir = os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', 'matrices', 'json'))

    if not os.path.exists(default_matrix_dir):
        raise IOError('The default matrix directory does not exists at %s!' % default_matrix_dir)

    return default_matrix_dir


def load_default_environments():
    default_envs = {
        'environments': read_all_env_json(os.path.join(matrix_dir(), 'all_matrices.json')),
        'log_pam1': read_env_json((os.path.join(matrix_dir(), 'logPAM1.json')))
    }
    return default_envs


def normalize_sequence(s):
    """
    Subtracts the ASCII value of 'A' from every character in the given string. This is necessary because the C core
    is working on these kind of inputs.

    :param s: the string which we would like to transform
     the underscore characters will not be transformed
    :return: the transformed string
    """
    ret = ''
    reg = re.compile('^[A-Z_]*$')

    if not reg.match(s):
        raise ValueError("Could not normalize '%s', because it contains invalid characters." % s)

    for c in s:
        if c != '_':
            ret += chr(ord(c) - ord('A'))
        else:
            #print "The '_' character has not been normalized!"
            ret += c

    return ret


def scale_to_byte(val, factor):
    """
    Scales the given double value to a byte value. This is needed to create the matrix and the gapOpen/Ext costs
    from the double versions, for the byte alignment function

    :param val: the double value which we want to scale up
    :param factor: the factor used for the scaling
    :return: a scaled integer from the -128..127 interval
    """
    ret = val * factor

    #for an upper bound estimation overflow is a problem, but underflow is not
    if ret > 127:
        raise OverflowError("Scaling overflow in scale_to_byte factor = %f, doubleValue = %f, result = %f" %
                (factor, val, ret))

    if ret < -128:
        return -128

    return int(math.ceil(ret))


def scale_to_short(val, factor):
    """
    Scales the given double value to a short value. This is needed to create the matrix and the gapOpen/Ext costs
    from the double versions, for the short alignment function

    :param val: the double value which we want to scale up
    :param factor: the factor used for the scaling
    :return: a scaled integer from the -32768..32767 interval
    """
    ret = val * factor

    #for an upper bound estimation overflow is a problem, but underflow is not
    if ret > 32767:
        raise OverflowError("Scaling overflow in scaleToShort factor = %f, doubleValue = %f, result = %f" %
                        (factor, val, ret))

    if ret < -32768:
        return -32768


    return int(math.ceil(ret))


def scale_back(val, factor):
    """
    Scales a double value back to a double by using the given factor.
    This is needed because the C core that returns a scaled result due to he scaled matrix and gap costs. We should use
    this function on the value returned by the short and byte alignment function written in C.

    :param val: the double  we want to scale down
    :param factor: the factor used for the scaling
    :return: a scaled double
    """
    return val / factor


def create_environment(gap_open, gap_ext, pam_distance, scores, column_order, threshold=85.0, **kwargs):
    """
    Creates an environment from the given parameters.
    :param gap_open: gap opening cost
    :param gap_ext: gap extension cost
    :param pam_distance: pam distance
    :param scores: distance matrix
    :param column_order: column order of the distance matrix
    :param kwargs:
    :return: an AlignmentEnvironment
    """

    reg = re.compile('^[A-Z]*$')
    column_order = ''.join(column_order)

    if not reg.match(column_order):
        raise ValueError("Could not create environment with columns '%s', because it contains invalid characters." %
                         column_order)

    if len(scores) != len(column_order):
        raise ValueError('The dimension of the matrix is not consistent with the column order')

    #TODO check whether gap_open <= gap_ext?


    env = AlignmentEnvironment()
    env.threshold = threshold
    env.gap_open = gap_open
    env.gap_ext = gap_ext
    env.pam = pam_distance
    compact_matrix = scores

    #convert the compact matrix into C compatible one by extending it to a 26x26 one
    extended_matrix = [[0 for x in xrange(26)] for x in xrange(26)]
    for i in range(0, len(column_order)):
        if len(compact_matrix[i]) != len(column_order):
            raise ValueError('The dimension of the matrix is not consistent with the column order')
        for j in range(0, len(column_order)):
            extended_matrix[ord(column_order[i]) - ord('A')][ord(column_order[j]) - ord('A')] = compact_matrix[i][j]

    env.float64_matrix = np.array(extended_matrix, dtype=np.float64)

    env.create_scaled_matrices()

    return env


def read_env_json(json_data):
    """
    This function reads an AlignmentEnvironment from a JSON object or a file that contains the JSON data
    :param json_data: the JSON object from which we want to read the environment or a JSON file
    :return: the environment
    """
    if isinstance(json_data, str):
        with open(json_data) as f:
            json_data = json.load(f)

    return create_environment(**json_data)


def read_all_env_json(file_loc):
    """
    Reads all of the alignment environments from the given json file.

    :param file_loc: the location where the matrices are stored
    :return: a list of AlignmentEnvironments
    """
    with open(file_loc) as json_file:
        json_data = json.load(json_file)

    ret = []

    for matrix_json in json_data['matrices']:
        ret.append(read_env_json(matrix_json))

    return ret


def align_short(s1, s2, env):
    """
    Aligns two sequences by using the matrix and gap costs defined in the file/AlignmentEnvironment.
    This is not an efficient way to align
    sequences and should only be used for testing purposes.

    :param s1: first string of the alignment
    :param s2: second string of the alignment
    :param env: file that contains the matrix in JSON format or an AlignmentEnvironment
    :return: the short estimation of the score
    """
    p = AlignmentProfile()
    if isinstance(env, str):
        env = read_env_json(env)
    p.create_profiles(s1, env)

    return p.align_short(s2, env)


def align_byte(s1, s2, env):
    """
    Aligns two sequences by using the matrix and gap costs defined in the file/AlignmentEnvironment.
    This is not an efficient way to align
    sequences and should only be used for testing purposes.

    :param s1: first string of the alignment
    :param s2: second string of the alignment
    :param env: file that contains the matrix in JSON format or an AlignmentEnvironment
    :return: the byte estimation of the score
    """
    p = AlignmentProfile()
    if isinstance(env, str):
        env = read_env_json(env)
    p.create_profiles(s1, env)

    return p.align_byte(s2, env)


cpdef double c_align_scalar_normalized_reference_local(np.ndarray[np.double_t,ndim=2] matrix, const char *s1,
                                                                int ls1, const char *s2, int ls2, double gap_open,
                                                                double gap_ext, double threshold):
    """
    This is a simple wrapper for the scalar alignment C function.
    :param matrix: the 26x26 double matrix
    :param s1: the first string
    :param ls1: the length of the first string
    :param s2: the second string
    :param ls2: the length of the second string
    :param gap_open: the gap opening cost
    :param gap_ext: the gap extension cost
    :param threshold: the threshold used for the calculation (might be ignored)
    :return: the exact scalar score
    """
    return cpyopa.c_align_scalar_reference_local(<double*> matrix.data, s1, ls1, s2, ls2,
                                                gap_open, gap_ext, threshold)


cpdef c_align_double_normalized_global(np.ndarray[np.double_t,ndim=2] matrix, const char *s1, int ls1, const char *s2,
                                       int ls2, double gap_open, double gap_ext):
    """
    A vectorized double precision GLOBAL alignment implementation.
    :param matrix: the 26x26 double matrix
    :param s1: the first string
    :param ls1: the length of the first string
    :param s2: the second string
    :param ls2: the length of the second string
    :param gap_open: the gap opening cost
    :param gap_ext: the gap extension cost
    :return: the score and the max1, max2 ranges
    """

    ret = []

    res = cpyopa.c_align_double_global(<double*> matrix.data, s1, ls1, s2, ls2,
                                      gap_open, gap_ext)

    ret.append(res)
    ret.append(ls1 - 1)
    ret.append(ls2 - 1)

    return ret


cpdef align_strings(s1, s2, env, is_global=False, provided_alignment=None):
    """
    Does the concrete string alignment for two given strings.
    :param s1: first sequence
    :param s2: second sequence
    :param env: AlignmentEnvironment to be used
    :param is_global:
    :param provided_alignment: If given, we can save up the calculation of the ranges and the scores.
     The provided alignment must contain the FULL ranges and the score.
    :return: the two aligned strings
    """
    if provided_alignment is None:
        provided_alignment = align_double(s1, s2, env, False, is_global, True)
    elif len(provided_alignment) != 5:
        raise ValueError('The provided alignment is invalid.'
                        ' It should contain the score, and the ranges for both sequences.')

    s1 = s1.s_norm
    s2 = s2.s_norm

    if not is_global:
        s1 = s1[provided_alignment[3]:provided_alignment[1] + 1]
        s2 = s2[provided_alignment[4]:provided_alignment[2] + 1]

    #TODO remove these char arrays
    cdef char o1[100010]
    cdef char o2[100010]

    aligned_s1 = ''
    aligned_s2 = ''

    cdef np.ndarray[np.double_t,ndim=2] matrix = env.float64_matrix

    max_len = cpyopa.c_align_strings(<double*> matrix.data, s1.encode('utf-8'), len(s1),
                                    s2.encode('utf-8'), len(s2), provided_alignment[0], o1, o2, 0.5e-4, env.gap_open, env.gap_ext)

    o1_p = [0] * max_len
    o2_p = [0] * max_len
    for i in range(max_len):
        o1_p[i] = o1[i]
        o2_p[i] = o2[i]

    aligned_s1 = Sequence(o1_p, True)
    aligned_s2 = Sequence(o2_p, True)

    return [aligned_s1, aligned_s2]

def align_double(s1, s2, env, stop_at_threshold=False, is_global=False, calculate_ranges=False):
    """
    A vectorized implementation of the double precision alignment algorithm.
    :param s1:
    :param s2:
    :param env:
    :param stop_at_threshold: if True, we terminate on reaching the threshold and return the score just over it
    :param is_global:
    :param calculate_ranges: if True, we calculate the full ranges for the alignment
    :return: an array of [score, max1, max2] if calculate ranges is false and  [score, max1, max2, min1, min2] if it's
    true
    """
    profile1 = AlignmentProfile()
    profile1.create_profile_double(s1, env.float64_matrix)

    if is_global:
        res = c_align_double_normalized_global(env.float64_matrix, s1.s_norm.encode('utf-8'), len(s1), s2.s_norm.encode('utf-8'), len(s2),
                             env.gap_open, env.gap_ext)
        if calculate_ranges :
            res.extend([0, 0])
    else:
        res = profile1.align_double(s2, env, stop_at_threshold)

        if calculate_ranges:
                s1r = Sequence(s1.s_norm[res[1]::-1], True)
                s2r = Sequence(s2.s_norm[res[2]::-1], True)

                profile1r = AlignmentProfile()
                profile1r.create_profile_double(s1r, env.float64_matrix)

                reversed = profile1r.align_double(s2r, env, stop_at_threshold)
                assert( abs(res[0]-reversed[0]) <= abs(res[0])*1e-10 )

                res.extend([res[1] - reversed[1], res[2] - reversed[2]])

    return res


def align_scalar_reference_local(s1, s2, env):
    """
    This is a simpler interface to the scalar alignment function written in C by using the AlignmentEnvironment python
    class. This is a reference implementation and not vectorized.
    :param s1: first string
    :param s2: second string
    :param env: the AlignmentEnvironment to be used
    :return: the exact scalar score
    """
    return c_align_scalar_normalized_reference_local(env.float64_matrix, s1.s_norm.encode('utf-8'), len(s1), s2.s_norm.encode('utf-8'), len(s2),
                                     env.gap_open, env.gap_ext, env.threshold)

cpdef generate_env(log_pam1_env, new_pam, threshold=85.0):
    """
    Generates a new AlignmentEnvironment from the given log_pam1 environment by using the new_pam pam distance and
    the given threshold.
    :param log_pam1_env: The log_pam1 environment
    :param new_pam: The new pam number to be used
    :param threshold: The new threshold
    :return: a freshly generated AlignmentEnvironment
    """
    env = AlignmentEnvironment()
    env.threshold = threshold
    env.pam = new_pam

    cdef np.ndarray[np.double_t, ndim=2, mode="c"] logpam = log_pam1_env.float64_matrix
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] generated_matrix = env.float64_matrix
    cpyopa.CreateOrigDayMatrix(<double*>logpam.data, new_pam, <double*> generated_matrix.data)

    env.gap_ext = -1.3961
    env.gap_open = -37.64 + 7.434 * log10(new_pam)

    env.create_scaled_matrices()

    return env


def generate_all_env(log_pam1_env, env_num=1266, starting_pam=0.049449734348559203348, threshold=85.0):
    """
    Generates a list of environments starting from the given pam distance, by using the formula of
    starting_pam = min((1 + 1/45.0) * starting_pam, starting_pam + 1) for further environments.
    :param log_pam1_env: The log_pam1 environment
    :param env_num: number of environments to be generated
    :param starting_pam: the pam number of the first environment
    :param threshold: threshold
    :return:
    """
    envs = []

    for i in range(env_num):
        envs.append(generate_env(log_pam1_env, starting_pam, threshold))
        starting_pam = min((1 + 1/45.0) * starting_pam, starting_pam + 1)

    return envs


class Sequence:
    def __init__(self, s, is_normalized=False):
    
        if isinstance(s, str):
            if not is_normalized:
                self.s_norm = normalize_sequence(s)
            else:
                self.s_norm = s
        elif isinstance(s, list):
            self.s_norm = ''.join(list(map(chr, s)))
            if not is_normalized:
                self.s_norm = normalize_sequence(self.s_norm)
        else:
            raise ValueError('Cannot construct Sequence from the given type!')

    def __str__(self):
        return self.convert_readable()

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.s_norm == other.s_norm
        elif isinstance(other, str):
            return self.s_norm == normalize_sequence(other)
        else:
            return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)

        if result is NotImplemented:
            return result

        return not result

    def __len__(self):
        return len(self.s_norm)

    def convert_readable(self):
        readable = ''

        for c in self.s_norm:
            if c == '_':
                readable += '_'
            else:
                readable += chr(ord('A') + ord(c))

        return readable


cdef class AlignmentProfile:
    """
    This class is used to facilitate multiple alignment by using the same query and matrix with the C core SSE byte and
    short alignments. From the given query and matrix we can generate a so-called "profile" in C, which makes it
    possible to efficiently align multiple consecutive S sequences to the same query. We can either create a short or a
    byte profile respectively for the short and byte C core alignment.
    """

    #to store and free the profile generated from C
    cdef cpyopa.ProfileByte* _c_profileByte
    cdef cpyopa.ProfileShort* _c_profileShort
    cdef cpyopa.ProfileDouble* _c_profileDouble


    def __cinit__(self):
        self._c_profileByte = NULL
        self._c_profileShort = NULL


    def __dealloc__(self):
        if self._c_profileByte is not NULL:
            cpyopa.c_free_profile_byte_sse_local(self._c_profileByte)

        if self._c_profileShort is not NULL:
            cpyopa.c_free_profile_short_sse_local(self._c_profileShort)

        if self._c_profileDouble is not NULL:
            cpyopa.free_profile_double_sse(self._c_profileDouble)


    def create_profiles(self, query, env):
        """
        Creates the byte and the short profile by using the environment's short and byte matrices
        :param query: the query from which we want to create the profiles
        :param env: the environment that contains the scaled byte and short matrices
        """
        self.create_profile_byte(query, env.int8_matrix)
        self.create_profile_short(query, env.int16_matrix)
        self.create_profile_double(query, env.float64_matrix)


    cpdef create_profile_byte(self, query, np.ndarray[np.int8_t,ndim=2] matrix):
        """
        Creates a byte profile from the given query and matrix. The matrix
        must be scaled to the byte version.

        :param query: the query sequence which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to bytes)
        """

        if matrix.shape[0] != 26 or  matrix.shape[1] != 26:
            raise ValueError("Invalid matrix shape, the matrix must be 26x26.")

        query = query.s_norm

        if self._c_profileByte is not NULL:
            cpyopa.c_free_profile_byte_sse_local(self._c_profileByte)

        self._c_profileByte = cpyopa.c_create_profile_byte_sse_local(query.encode('utf-8'), len(query), <signed char*> matrix.data)


    cpdef create_profile_short(self, query, np.ndarray[np.int16_t,ndim=2] matrix):
        """
        Creates a short profile from the given query and matrix. The matrix
        must be scaled to the short version.

        :param query: the query sequence which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to shorts)
        """

        if matrix.shape[0] != 26 or matrix.shape[1] != 26:
            raise ValueError("Invalid matrix shape, the matrix must be 26x26.")

        query = query.s_norm

        if self._c_profileShort is not NULL:
            cpyopa.c_free_profile_short_sse_local(self._c_profileShort)

        self._c_profileShort = cpyopa.c_create_profile_short_sse_local(query.encode('utf-8'), len(query),
                                                                      <signed short*> matrix.data)


    cpdef create_profile_double(self, query, np.ndarray[np.double_t, ndim=2] matrix):
        """
        Creates a double profile from the given query and matrix.

        :param query: the query sequence which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to shorts)
        """

        if matrix.shape[0] != 26 or matrix.shape[1] != 26:
            raise ValueError("Invalid matrix shape, the matrix must be 26x26.")

        query = query.s_norm

        if self._c_profileDouble is not NULL:
            cpyopa.free_profile_double_sse(self._c_profileDouble)

        self._c_profileDouble = cpyopa.createProfileDoubleSSE(query.encode('utf-8'), len(query), <double*> matrix.data)


    cpdef align_byte(self, s2, env):
        """
        Aligns the sequence s2 with the gap costs defined in the given AlignmentEnvironment to the profile by using byte
        SSE alignment.

        :param s2: the sequence which we would like to align to the profile
        :param env: the AlignmentEnvironment that contains the gap opening/extending costs and the threshold
        :return: the result of the byte alignment
        """

        s2 = s2.s_norm

        ret = cpyopa.c_align_profile_byte_sse_local(self._c_profileByte, s2.encode('utf-8'), len(s2),
                                                   env.int8_gap_open, env.int8_gap_ext, env.threshold)
        #do not scale back DBL_MAX
        if ret >= sys.float_info.max:
            return ret

        return scale_back(ret, env.byte_factor())


    cpdef align_short(self, s2, env):
        """
        Aligns the sequence s2 with the gap costs defined in the given AlignmentEnvironment to the profile
        by using short SSE alignment.

        :param env: the AlignmentEnvironment that contains the gap opening/extending costs and the threshold
        :return: the result of the short alignment
        """

        s2 = s2.s_norm

        ret = cpyopa.c_align_profile_short_sse_local(self._c_profileShort,
                                                    s2.encode('utf-8'), len(s2), env.int16_gap_open, env.int16_gap_ext, env.threshold)
        #do not scale back DBL_MAX
        if ret >= sys.float_info.max:
            return ret

        return scale_back(ret, env.short_factor())


    cpdef align_double(self, s2, env, stop_at_threshold=False):
        """
        Aligns the sequence s2 with the gap costs defined in the given AlignmentEnvironment
        to the profile by using double SSE alignment.

        :param s2: the sequence which we would like to align to the profile
        :param env: the AlignmentEnvironment that contains the gap opening/extending costs and the threshold
        :return: the result of the double alignment
        """
        cdef int max1[1]
        cdef int max2[1]

        max1[0] = 0
        max2[0] = 0

        ret = []
        s2 = s2.s_norm

        if not stop_at_threshold:
            res = cpyopa.align_double_local(self._c_profileDouble, s2.encode('utf-8'), len(s2),
                                           env.gap_open, env.gap_ext, sys.float_info.max, max1, max2)
        else:
            res = cpyopa.align_double_local(self._c_profileDouble, s2.encode('utf-8'), len(s2),
                                           env.gap_open, env.gap_ext, env.threshold, max1, max2)

        ret.append(res)
        ret.append(max1[0] - 1)
        ret.append(max2[0] - 1)

        return ret


class AlignmentEnvironment:
    """
    This class stores all the information that is necessary to conduct an alignment (including the matrix and the gap
    costs). It does the necessary transformations such as scaling up the matrix and gap costs and extending the matrix
    to a 26x26 one, which is needed for scalar and byte alignments. Although you can manually create the necessary
    matrices and gap costs it is highly recommended to use this class for alignments.
    """
    def __init__(self):
        #the PAM distance associated with the matrices
        self.pam = 0.0

        self.threshold = 85.0

        #This gapOpen/Ext and matrix is used by the scalar alignment function. These are directly passed to the C core
        self.gap_open = 0.0
        self.gap_ext = 0.0
        self.float64_matrix = np.zeros(shape=(26, 26), dtype=np.float64)

        #This gapOpen/Ext and matrix is used by the short alignment function. These are directly passed to the C core
        #and calculated from the scalar matrix and gap costs (by a simple scaling)
        self.int16_matrix = np.zeros(shape=(26, 26), dtype=np.int16)
        self.int16_gap_open = 0.0
        self.int16_gap_ext = 0.0

        #This gapOpen/Ext and matrix is used by the byte alignment function. These are directly passed to the C core
        #and calculated from the scalar matrix and gap costs (by a simple scaling)
        self.int8_matrix = np.zeros(shape=(26, 26), dtype=np.int8)
        self.int8_gap_open = 0.0
        self.int8_gap_ext = 0.0


    def create_scaled_matrices(self):
        """
        Creates the int8 and int16 matrix and gap costs from the double matrix by using a simple scaling
        :return:nothing
        """

        factor = self.byte_factor()
        self.int8_matrix = np.vectorize(lambda x: scale_to_byte(x, factor))(self.float64_matrix ).astype(np.int8)
        self.int8_gap_open = scale_to_byte(self.gap_open, factor)
        self.int8_gap_ext = scale_to_byte(self.gap_ext, factor)

        factor = self.short_factor()
        self.int16_matrix = np.vectorize(lambda x: scale_to_short(x, factor))(self.float64_matrix ).astype(np.int16)
        self.int16_gap_open = scale_to_short(self.gap_open, factor)
        self.int16_gap_ext = scale_to_short(self.gap_ext, factor)

    #calculates the byte scaling factor
    def byte_factor(self):

        #This is copied from the C code and I have no idea why we use this at the byte version but not at the short
        abs_min = abs(np.amin(self.float64_matrix))
        return 255.0/(self.threshold + abs_min)


    def short_factor(self):
        return 65535.0 / self.threshold



cdef class MutipleAlEnv:
    """
    This class stores a list of DayHoff matrices generated from C. It can be generated from a list of
     AlignmentEnvironments and can be used to call the EstimatePam function, which is implemented in C.
    """

    cdef cpyopa.DayMatrix* _c_dayMatrices
    cdef int dms_len
    cdef log_pam1

    def __init__(self, envs, log_pam1):
        self.create_day_matrices(envs)
        self.log_pam1 = log_pam1

    def __cinit__(self):
        self._c_dayMatrices = NULL

    def __dealloc__(self):
        self.free()

    def create_day_matrices(self, envs):
        self.free()

        self.dms_len = len(envs)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] gap_open = np.empty(shape=(self.dms_len + 1))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] gap_ext = np.empty(shape=(self.dms_len + 1))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] pam_dist = np.empty(shape=(self.dms_len + 1))
        cdef np.ndarray[np.uint64_t, ndim=1, mode="c"] matrices = np.empty(shape=(self.dms_len + 1), dtype=np.uint64)

        for i in range(1, self.dms_len + 1):
            gap_open[i] = envs[i-1].gap_open
            gap_ext[i] = envs[i-1].gap_ext
            pam_dist[i] = envs[i-1].pam
            matrices[i] = envs[i-1].float64_matrix.ctypes.data

        self._c_dayMatrices = cpyopa.createDayMatrices(<double*> gap_open.data, <double*> gap_ext.data,
                                                      <double*> pam_dist.data, <long long*> matrices.data, self.dms_len)

    def estimate_pam(self, s1, s2):
        if len(s1) != len(s2):
            raise ValueError('The length of s1 and s2 must be equal!')

        cdef double result[3]
        cdef np.ndarray[np.double_t, ndim=2, mode="c"] logpam = self.log_pam1.float64_matrix

        cpyopa.EstimatePam(s1.s_norm.encode('utf-8'), s2.s_norm.encode('utf-8'), len(s1),
                          self._c_dayMatrices, self.dms_len, <double*> logpam.data, result)
        ret = [0.0] * 3
        ret[0] = result[0]
        ret[1] = result[1]
        ret[2] = result[2]

        return ret

    def free(self):
        if self._c_dayMatrices is not NULL:
            cpyopa.freeDayMatrices(self._c_dayMatrices, self.dms_len)

