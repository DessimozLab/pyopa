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


def readAlignmentEnvironments(loc):
    ret = []
    listDir = sorted(os.listdir(loc), key=lambda x: int(os.path.splitext(os.path.basename(x))[0]))

    for fn in listDir:
        f = os.path.join(loc, fn)
        if os.path.isfile(f) and fnmatch.fnmatch(f, "*.dat"):
            tmp = AlignmentEnvironment()
            tmp.readFromFile(f)
            ret.append(tmp)

    return ret


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


    cpdef createProfileByte(self, query, np.ndarray[np.int8_t,ndim=2] matrix):
        """
        Creates a byte profile from the given query and matrix. The query string must not be normalized but the matrix
        must be scaled to the byte version.

        :param query: the query string which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to bytes)
        """

        if matrix.shape[0] != 26 or  matrix.shape[1] != 26:
            raise Exception("Invalid matrix shape, the matrix must be 26x26.")

        if self._c_profileByte is not NULL:
            cython_swps3.python_freeProfileByteSSE(self._c_profileByte)

        q = normalizeString(query)
        self._c_profileByte = cython_swps3.python_createByteProfileSSE(q, len(query), <signed char*> matrix.data)


    cpdef createProfileShort(self, query, np.ndarray[np.int16_t,ndim=2] matrix):
        """
        Creates a short profile from the given query and matrix. The query string must not be normalized but the matrix
        must be scaled to the short version.

        :param query: the query string which we want to use for the profile
        :param matrix: the matrix which we want to use for the profile (must be scaled to shorts)
        """

        if matrix.shape[0] != 26 or matrix.shape[1] != 26:
            raise Exception("Invalid matrix shape, the matrix must be 26x26.")

        if self._c_profileShort is not NULL:
            cython_swps3.python_freeProfileShortSSE(self._c_profileShort)

        q = normalizeString(query)
        self._c_profileShort = cython_swps3.python_createShortProfileSSE(q, len(query), <signed short*> matrix.data)


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
        #TODO replace this constant with a better solution
        if ret >= 1.7976931348623157e+308:
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
        #TODO replace this constant with a better solution
        if ret >= 1.7976931348623157e+308:
            return ret

        return scaleBack(ret, env.shortFactor())

    def alignByte(self, s2, env):
        """
        Aligns s2 to the profile with byte alignment. This function automatically normalizes the s2 input string. This
        function is not efficient because we have to transform the input string on every call.

        :param s2: the string to which we would like to align the profile (all characters must be in range of A-Z)
        :param env: the AlignmentEnvironment that contains the gap costs and the threshold
        :return: the result of the byte alignment
        """
        return self.alignByteNormalized(normalizeString(s2), env)

    def alignShort(self, s2, env):
        """
        Aligns s2 to the profile with short alignment. This function automatically normalizes the s2 input string. This
        function is not efficient because we have to transform the input string on every call.

        :param s2: the string to which we would like to align the profile (all characters must be in range of A-Z)
        :param env: the AlignmentEnvironment that contains the gap costs and the threshold
        :return: the result of the short alignment
        """
        return self.alignShortNormalized(normalizeString(s2), env)


class AlignmentEnvironment:
    def __init__(self):
        self.columns = "ARNDCQEGHILKMFPSTWYVX"
        self.gapOpen = 0.0
        self.gapExt = 0.0
        self.threshold = 85.0
        self.float64_Matrix = np.ndarray(shape=(26, 26), dtype=np.float64)
        self.matrixId = ""

        self.int16_Matrix = np.ndarray(shape=(26, 26), dtype=np.int16)
        self.int16_gapOpen = 0.0
        self.int16_gapExt = 0.0

        self.int8_Matrix = np.ndarray(shape=(26, 26), dtype=np.int8)
        self.int8_gapOpen = 0.0
        self.int8_gapExt = 0.0


    def readFromFile(self, file):
        self.matrixId = os.path.splitext(os.path.basename(file))[0]

        with open(file) as data_file:    
            data = json.load(data_file)

        self.gapOpen = data["GapOpen"]
        self.gapExt = data["GapExt"]
        compactMatrix = data["Scores"]

        #convert the compact matrix into C compatible one
        extendedMatrix = [[0 for x in xrange(26)] for x in xrange(26)]
        for i in range(0, len(self.columns)):        
            for j in range(0, len(self.columns)):
                extendedMatrix[ord(self.columns[i]) - ord('A')][ord(self.columns[j]) - ord('A')] = compactMatrix[i][j]


        #align NumPy array pointer to 16 byte
        self.float64_Matrix = np.array(extendedMatrix, dtype=np.float64)

        #create the int8 and int16 matrix from the double matrix by using a simple scaling
        factor = self.byteFactor()
        self.int8_Matrix = np.vectorize(lambda x: scaleToByte(x, factor))(self.float64_Matrix ).astype(np.int8)
        self.int8_gapOpen =  scaleToByte(self.gapOpen, factor)
        self.int8_gapExt =  scaleToByte(self.gapExt, factor)

        factor = self.shortFactor()
        self.int16_Matrix = np.vectorize(lambda x: scaleToShort(x, factor))(self.float64_Matrix ).astype(np.int16)
        self.int16_gapOpen =  scaleToShort(self.gapOpen, factor)
        self.int16_gapExt =  scaleToShort(self.gapExt, factor)


    def byteFactor(self):
        absMin = abs(np.amin(self.float64_Matrix))
        return 255.0/(self.threshold + absMin)


    def shortFactor(self):
        return 65535.0 / self.threshold
