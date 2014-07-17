cimport cython_swps3

import numpy as np
cimport numpy as np
import json
import sys
import math
import os
import fnmatch


def normalizeString(s):
    ret = ""
    for c in s:
        ret+=chr(ord(c) - ord('A'))

    return ret


def scaleToByte(dVal, factor):
    ret = dVal * factor

    if ret > 127 or ret < -128:
        print "Error: scaling overflow in scaleToByte factor = %f, doubleValue = %f" % (factor, dVal)

    return int(math.ceil(ret))


def scaleBackFromByte(bVal, factor):
    return bVal / factor


def scaleBackFromShort(sVal, factor):
    return sVal / factor


def scaleToShort(dVal, factor):
    ret = dVal * factor

    if ret > 32767 or ret < -32768:
        print "Error: scaling overflow in scaleToShort factor = %f, doubleValue = %f" % (factor, dVal)

    return int(math.ceil(ret))


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


cdef class Profile:
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


    cpdef createProfileByte(self, query, np.ndarray[np.int8_t,ndim=2] matrix):
        if self._c_profileByte is not NULL:
            cython_swps3.python_freeProfileByteSSE(self._c_profileByte)

        self._c_profileByte = cython_swps3.python_createByteProfileSSE(query, len(query), <signed char*> matrix.data)


    cpdef createProfileShort(self, query, np.ndarray[np.int16_t,ndim=2] matrix):
        if self._c_profileShort is not NULL:
            cython_swps3.python_freeProfileShortSSE(self._c_profileShort)

        self._c_profileShort = cython_swps3.python_createShortProfileSSE(query, len(query), <signed short*> matrix.data)


    cpdef alignByteProfile(self, s2, env):
        ret = cython_swps3.python_alignByteProfileSSE(<ProfileByte*> self._c_profileByte, s2, len(s2), env.int8_gapOpen, env.int8_gapExt, env.threshold)
        #do not scale back DBL_MAX
        if ret >= 1.7976931348623157e+308:
            return ret
        return scaleBackFromByte(ret, env.byteFactor())


    cpdef alignShortProfile(self, s2, env):
        ret = cython_swps3.python_alignShortProfileSSE(<ProfileShort*> self._c_profileShort, s2, len(s2), env.int16_gapOpen, env.int16_gapExt, env.threshold)
        #do not scale back DBL_MAX
        if ret >= 1.7976931348623157e+308:
            return ret
        return scaleBackFromShort(ret, env.shortFactor())



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
