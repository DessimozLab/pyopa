cimport cython_swps3

import numpy as np
cimport numpy as np
import json
import sys
import math
import os
import fnmatch

#This function expects "normalized" byte arrays as input parameters, which means that 'A' should be normalized to 0
# 'B' to 1, 'C' to 2 etc... and by 0, 1, 2 NOT character '0', '1', '2' is meant
cpdef double swps3_alignScalar(np.ndarray[np.double_t,ndim=2] matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    return cython_swps3.python_alignScalar(<double*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


#This function expects "normalized" byte arrays as input parameters, which means that 'A' should be normalized to 0
# 'B' to 1, 'C' to 2 etc... and by 0, 1, 2 NOT character '0', '1', '2' is meant
#
# In addition, the gapOpen, gapExt, and the Matrix should also be scaled by a factor of SIGNED_CHAR_MAX / THRESHOLD
# The returned double is NOT scaled back by the factor, so a scaleback is necessary
cpdef double swps3_alignByte(np.ndarray[np.int8_t,ndim=2] matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    return cython_swps3.python_alignByte(<signed char*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


#This function expects "normalized" byte arrays as input parameters, which means that 'A' should be normalized to 0
# 'B' to 1, 'C' to 2 etc... and by 0, 1, 2 NOT character '0', '1', '2' is meant
#
# In addition, the gapOpen, gapExt, and the Matrix should also be scaled by a factor of SIGNED_SHORT_MAX / THRESHOLD
# The returned double is NOT scaled back by the factor, so a scaleback is necessary
cpdef double swps3_alignShort(np.ndarray[np.int8_t,ndim=2] matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    #TODO this might not work correctly because C only contains byte matrix implementation???
    #TODO the matrix should not be casted to int8
    return cython_swps3.python_alignShort(<signed char*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


def alignScalar(s1, s2, env):
    return swps3_alignScalar(env.float64_Matrix, s1, len(s1), s2, len(s2), env.gapOpen, env.gapExt, env.threshold)


def normalizedAlignScalar(s1, s2, env):
    return alignScalar(normalizeString(s1), normalizeString(s2), env)


def alignByte(s1, s2, env):
    ret = swps3_alignByte(env.int8_Matrix, s1, len(s1), s2, len(s2), env.int8_gapOpen, env.int8_gapExt, env.threshold)
    return scaleBackFromByte(ret, env.byteFactor())


def normalizedAlignByte(s1, s2, env):
    return alignByte(normalizeString(s1), normalizeString(s2), env)


#TODO should not transform into int8
def alignShort(s1, s2, env):
    tmp = env.int16_Matrix.astype(np.int8)
    ret = swps3_alignShort(tmp, s1, len(s1), s2, len(s2), env.int16_gapOpen, env.int16_gapExt, env.threshold)
    return scaleBackFromShort(ret, env.shortFactor())

def normalizedAlignShort(s1, s2, env):
    return alignShort(normalizeString(s1), normalizeString(s2), env)


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


#TODO this is linear for TESTING only
def binaryAlign(environments, s1, s2):
    minId = 0
    maxId = len(environments) - 1

    maxScore = alignScalar(s1, s2, environments[0]);

    for env in environments:
        tmp = alignScalar(s1, s2, env)

        if tmp > maxScore:
            maxScore = tmp

    return maxScore


def binaryAlignNormalized(environments, s1, s2):
    return binaryAlign(environments, normalizeString(s1), normalizeString(s2))


#TODO check the correctness of this!
def alignMatrix16Byte(matrix):
    if len(matrix.shape) != 2:
        print "The dimension of the matrix must be 2!"
        return matrix

    n = matrix.shape[0]
    m = matrix.shape[1]

    dtype = matrix.dtype
    nbytes = n * m * dtype.itemsize
    buf = np.empty(nbytes + 16, dtype=np.uint8)
    start_index = -buf.ctypes.data % 16
    a = buf[start_index:start_index + nbytes].view(dtype).reshape(n, m)

    np.copyto(a, matrix)
    return a


class AlignmentEnvironment:
    def __init__(self):
        self.columns = "ARNDCQEGHILKMFPSTWYV"
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


        #TODO verify the correctness of the 16 byte alignment
        #align NumPy array pointer to 16 byte
        self.float64_Matrix = alignMatrix16Byte(np.array(extendedMatrix, dtype=np.float64))

        #create the int8 and int16 matrix from the double matrix by using a simple scaling
        factor = self.byteFactor()
        self.int8_Matrix = alignMatrix16Byte(np.vectorize(lambda x: scaleToByte(x, factor))(self.float64_Matrix ).astype(np.int8))
        self.int8_gapOpen =  scaleToByte(self.gapOpen, factor)
        self.int8_gapExt =  scaleToByte(self.gapExt, factor)

        factor = self.shortFactor()
        self.int16_Matrix = alignMatrix16Byte(np.vectorize(lambda x: scaleToShort(x, factor))(self.float64_Matrix ).astype(np.int16))
        self.int16_gapOpen =  scaleToShort(self.gapOpen, factor)
        self.int16_gapExt =  scaleToShort(self.gapExt, factor)


    def byteFactor(self):
        absMin = abs(np.amin(self.float64_Matrix))
        return 255.0/(self.threshold + absMin)


    def shortFactor(self):
        return 65535.0 / self.threshold
