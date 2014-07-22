cimport cython_swps3

import numpy as np
cimport numpy as np
import json
import sys
import math
import os
import fnmatch
import re


def normalizeString(s):
    """
    Subtracts the ASCII value of 'A' from every character in the given string. This is necessary because the C core
    is working on these kind of inputs. Normally you should not call this function since it is automatically called
    during the profile generation.

    :param s: the string which we would like to transform
    :return: the transformed string
    """
    ret = ""
    reg = re.compile('^[A-Z]*$')

    if not reg.match(s):
        raise Exception("Could not normalize '%s', because it contains invalid characters." % s)

    for c in s:
        ret+=chr(ord(c) - ord('A'))

    return ret


def scaleToByte(dVal, factor):
    """
    Scales the given double value to a byte value. This is needed to create the matrix and the gapOpen/Ext costs
    from the double versions, for the byte alignment function

    :param dVal: the double value which we want to scale up
    :param factor: the factor used for the scaling
    :return: a scaled integer from the -128..127 interval
    """
    ret = dVal * factor

    if ret > 127 or ret < -128:
        raise Exception("Scaling overflow in scaleToByte factor = %f, doubleValue = %f, result = %f" %
                        (factor, dVal, ret))

    return int(math.ceil(ret))


def scaleToShort(dVal, factor):
    """
    Scales the given double value to a short value. This is needed to create the matrix and the gapOpen/Ext costs
    from the double versions, for the short alignment function

    :param dVal: the double value which we want to scale up
    :param factor: the factor used for the scaling
    :return: a scaled integer from the -32768..32767 interval
    """
    ret = dVal * factor

    if ret > 32767 or ret < -32768:
        raise Exception("Scaling overflow in scaleToShort factor = %f, doubleValue = %f, result = %f" %
                        (factor, dVal, ret))

    return int(math.ceil(ret))


def scaleBack(val, factor):
    """
    Scales a double value back to a double by using the given factor.
    This is needed because the C core that returns a scaled result due to he scaled matrix and gap costs. We should use
    this function on the value returned by the short and byte alignment function written in C.

    :param val: the double  we want to scale down
    :param factor: the factor used for the scaling
    :return: a scaled double
    """
    return val / factor


def readEnvJson(file):
    """
    This function is reading an AlignmentEnvironment from a JSON formatted file.
    :param file: the file from which we want to read the environment
    :return: the environment
    """

    env = AlignmentEnvironment()

    with open(file) as data_file:
        data = json.load(data_file)

    env.gapOpen = data["GapOpen"]
    env.gapExt = data["GapExt"]
    env.pam = data["PamDistance"]
    compactMatrix = data["Scores"]

    #convert the compact matrix into C compatible one by extending it to a 26x26 one
    extendedMatrix = [[0 for x in xrange(26)] for x in xrange(26)]
    for i in range(0, len(env.columns)):
        for j in range(0, len(env.columns)):
            extendedMatrix[ord(env.columns[i]) - ord('A')][ord(env.columns[j]) - ord('A')] = compactMatrix[i][j]


    env.float64_Matrix = np.array(extendedMatrix, dtype=np.float64)

    #create the int8 and int16 matrix from the double matrix by using a simple scaling
    factor = env.byteFactor()
    env.int8_Matrix = np.vectorize(lambda x: scaleToByte(x, factor))(env.float64_Matrix ).astype(np.int8)
    env.int8_gapOpen =  scaleToByte(env.gapOpen, factor)
    env.int8_gapExt =  scaleToByte(env.gapExt, factor)

    factor = env.shortFactor()
    env.int16_Matrix = np.vectorize(lambda x: scaleToShort(x, factor))(env.float64_Matrix ).astype(np.int16)
    env.int16_gapOpen =  scaleToShort(env.gapOpen, factor)
    env.int16_gapExt =  scaleToShort(env.gapExt, factor)

    return env


def readAllEnvJson(loc, ext):
    """
    Reads all of the alignment environments from the given directory. The extension of the files must be .dat

    :param loc: the directory where the matrices are stored
    :return: a list of AlignmentEnvironments sorted by the file name
    """
    ret = []
    listDir = sorted(os.listdir(loc), key=lambda x: int(os.path.splitext(os.path.basename(x))[0]))

    for fn in listDir:
        f = os.path.join(loc, fn)
        if os.path.isfile(f) and fnmatch.fnmatch(f, "*." + ext):
            ret.append(readEnvJson(f))

    return ret


def alignShort(s1, s2, file):
    """
    Aligns two sequences by using the matrix and gap costs defined in the file. This is not an efficient way to align
    sequences and should only be used for testing purposes.

    :param s1: first string of the alignment
    :param s2: second string of the alignment
    :param file: file that contains the matrix in JSON format
    :return: the short estimation of the score
    """
    p = AlignmentProfile()
    e = readEnvJson(file)
    p.createProfiles(s1, e)

    return p.alignShort(s2, e)


def alignByte(s1, s2, file):
    """
    Aligns two sequences by using the matrix and gap costs defined in the file. This is not an efficient way to align
    sequences and should only be used for testing purposes.

    :param s1: first string of the alignment
    :param s2: second string of the alignment
    :param file: file that contains the matrix in JSON format
    :return: the byte estimation of the score
    """
    p = AlignmentProfile()
    e = readEnvJson(file)
    p.createProfiles(s1, e)

    return p.alignByte(s2, e)


cpdef double alignScalarNormalizedC(np.ndarray[np.double_t,ndim=2] matrix, const char *s1, int ls1,
                                    const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    """
    This is a simple wrapper for the scalar alignment C function.
    :param matrix: the 26x26 double matrix
    :param s1: the first string
    :param ls1: the length of the first string
    :param s2: the second string
    :param ls2: the length of the second string
    :param gapOpen: the gap opening cost
    :param gapExt: the gap extension cost
    :param threshold: the threshold used for the calculation (might be ignored)
    :return: the exact scalar score
    """
    return cython_swps3.python_alignScalar(<double*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


def alignScalarNormalized(s1, s2, env):
    """
    Similar to alignScalar but working on normalized strings so more efficient.
    :param s1: first string
    :param s2: second string
    :param env: the AlignmentEnvironment to be used
    :return: the exact scalar score
    """
    return alignScalarNormalizedC(env.float64_Matrix, s1, len(s1), s2, len(s2), env.gapOpen, env.gapExt, env.threshold)


def alignScalar(s1, s2, env):
    """
    This is a simpler interface to the scalar alignment function written in C by using the AlignmentEnvironment python
    class. This is less efficient than the alignScalarNormalized version
    :param s1: first string
    :param s2: second string
    :param env: the AlignmentEnvironment to be used
    :return: the exact scalar score
    """
    ns1 = normalizeString(s1)
    ns2 = normalizeString(s2)
    return alignScalarNormalized(ns1, ns2, env)


cdef class AlignmentProfile:
    """
    This class is used to facilitate multiple alignment by using the same query and matrix with the C core SSE byte and
    short alignments. From the given query and matrix we can generate a so-called "profile" in C, which makes it
    possible to efficiently align multiple consecutive S strings to the same query. We can either create a short or a
    byte profile respectively for the short and byte C core alignment.
    """

    #to store and free the profile generated from C
    cdef cython_swps3.ProfileByte* _c_profileByte
    cdef cython_swps3.ProfileShort* _c_profileShort


    def __cinit__(self):
        self._c_profileByte = NULL
        self._c_profileShort = NULL


    def __dealloc__(self):
        if self._c_profileByte is not NULL:
            cython_swps3.python_freeProfileByteSSE(self._c_profileByte)

        if self._c_profileShort is not NULL:
            cython_swps3.python_freeProfileShortSSE(self._c_profileShort)


    def createProfiles(self, query, env):
        """
        Creates the byte and the short profile by using the environment's short and byte matrices
        :param query: the query from which we want to create the profiles
        :param env: the environment that contains the scaled byte and short matrices
        """
        self.createProfileByte(query, env.int8_Matrix)
        self.createProfileShort(query, env.int16_Matrix)


    def createProfilesNormalized(self, query, env):
        """
        Creates the byte and the short profile by using the environment's short and byte matrices
        :param query: the query from which we want to create the profiles
        :param env: the environment that contains the scaled byte and short matrices
        """
        self.createProfileByteNormalized(query, env.int8_Matrix)
        self.createProfileShortNormalized(query, env.int16_Matrix)


    cpdef createProfileByteNormalized(self, query, np.ndarray[np.int8_t,ndim=2] matrix):
        """
        Creates a byte profile from the given query and matrix. The query string must be normalized but the matrix
        must be scaled to the byte version.

        :param query: the query string which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to bytes)
        """

        if matrix.shape[0] != 26 or  matrix.shape[1] != 26:
            raise Exception("Invalid matrix shape, the matrix must be 26x26.")

        if self._c_profileByte is not NULL:
            cython_swps3.python_freeProfileByteSSE(self._c_profileByte)

        self._c_profileByte = cython_swps3.python_createByteProfileSSE(query, len(query), <signed char*> matrix.data)


    def createProfileByte(self, query, np.ndarray[np.int8_t,ndim=2] matrix):
        """
        Creates a byte profile from the given query and matrix. The query string must not be normalized but the matrix
        must be scaled to the byte version.

        :param query: the query string which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to bytes)
        """

        q = normalizeString(query)
        self.createProfileByteNormalized(q, matrix)


    cpdef createProfileShortNormalized(self, query, np.ndarray[np.int16_t,ndim=2] matrix):
        """
        Creates a short profile from the given query and matrix. The query string must be normalized but the matrix
        must be scaled to the short version.

        :param query: the query string which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to shorts)
        """

        if matrix.shape[0] != 26 or matrix.shape[1] != 26:
            raise Exception("Invalid matrix shape, the matrix must be 26x26.")

        if self._c_profileShort is not NULL:
            cython_swps3.python_freeProfileShortSSE(self._c_profileShort)

        self._c_profileShort = cython_swps3.python_createShortProfileSSE(query, len(query), <signed short*> matrix.data)


    def createProfileShort(self, query, np.ndarray[np.int16_t,ndim=2] matrix):
        """
        Creates a short profile from the given query and matrix. The query string must not be normalized but the matrix
        must be scaled to the short version.

        :param query: the query string which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to shorts)
        """

        q = normalizeString(query)
        self.createProfileByteNormalized(q, matrix)


    cpdef alignByteNormalized(self, s2, env):
        """
        Aligns the string s2 with the gap costs defined in the given AlignmentEnvironment to the profile by using byte
        SSE alignment.

        :param s2: the string which we would like to align to the profile, this input string must be normalized
        :param env: the AlignmentEnvironment that contains the gap opening/extending costs and the threshold
        :return: the result of the byte alignment
        """
        ret = cython_swps3.python_alignByteProfileSSE(<ProfileByte*> self._c_profileByte, s2, len(s2),
                                                      env.int8_gapOpen, env.int8_gapExt, env.threshold)
        #do not scale back DBL_MAX
        if ret >= sys.float_info.max:
            return ret

        return scaleBack(ret, env.byteFactor())


    cpdef alignShortNormalized(self, s2, env):
        """
        Aligns the string s2 with the gap costs defined in the given AlignmentEnvironment to the profile by using short
        SSE alignment.

        :param s2: the string which we would like to align to the profile, this input string must be normalized
        :param env: the AlignmentEnvironment that contains the gap opening/extending costs and the threshold
        :return: the result of the short alignment
        """
        ret = cython_swps3.python_alignShortProfileSSE(<ProfileShort*> self._c_profileShort,
                                                       s2, len(s2), env.int16_gapOpen, env.int16_gapExt, env.threshold)
        #do not scale back DBL_MAX
        if ret >= sys.float_info.max:
            return ret

        return scaleBack(ret, env.shortFactor())

    def alignByte(self, s2, env):
        """
        Aligns s2 to the profile with byte alignment. This method automatically normalizes the s2 input string. This
        method is not efficient because we have to transform the input string on every call.

        :param s2: the string to which we would like to align the profile (all characters must be in range of A-Z)
        :param env: the AlignmentEnvironment that contains the gap costs and the threshold
        :return: the result of the byte alignment
        """
        return self.alignByteNormalized(normalizeString(s2), env)

    def alignShort(self, s2, env):
        """
        Aligns s2 to the profile with short alignment. This method automatically normalizes the s2 input string. This
        method is not efficient because we have to transform the input string on every call.

        :param s2: the string to which we would like to align the profile (all characters must be in range of A-Z)
        :param env: the AlignmentEnvironment that contains the gap costs and the threshold
        :return: the result of the short alignment
        """
        return self.alignShortNormalized(normalizeString(s2), env)


class AlignmentEnvironment:
    """
    This class stores all the information that is necessary to conduct an alignment (including the matrix and the gap
    costs). It does the necessary transformations such as scaling up the matrix and gap costs and extending the matrix
    to a 26x26 one, which is needed for scalar and byte alignments. Although you can manually create the necessary
    matrices and gap costs it is highly recommended to use this class for alignments.
    """
    def __init__(self):
        #the PAM distance associated with the matrices
        self.pam = 0

        #The columns of the compressed matrix. The compressed matrix to be extended to a 26x26 matrix since the profile
        #generation requires a 26x26 matrix
        self.columns = "ARNDCQEGHILKMFPSTWYV"
        self.threshold = 85.0

        #This gapOpen/Ext and matrix is used by the scalar alignment function. These are directly passed to the C core
        self.gapOpen = 0.0
        self.gapExt = 0.0
        self.float64_Matrix = np.ndarray(shape=(26, 26), dtype=np.float64)

        #This gapOpen/Ext and matrix is used by the short alignment function. These are directly passed to the C core
        #and calculated from the scalar matrix and gap costs (by a simple scaling)
        self.int16_Matrix = np.ndarray(shape=(26, 26), dtype=np.int16)
        self.int16_gapOpen = 0.0
        self.int16_gapExt = 0.0

        #This gapOpen/Ext and matrix is used by the byte alignment function. These are directly passed to the C core
        #and calculated from the scalar matrix and gap costs (by a simple scaling)
        self.int8_Matrix = np.ndarray(shape=(26, 26), dtype=np.int8)
        self.int8_gapOpen = 0.0
        self.int8_gapExt = 0.0


    #calculates the byte scaling factor
    def byteFactor(self):

        #This is copied from the C code and I have no idea why we use this at the byte version but not at the short
        absMin = abs(np.amin(self.float64_Matrix))
        return 255.0/(self.threshold + absMin)


    def shortFactor(self):
        return 65535.0 / self.threshold
