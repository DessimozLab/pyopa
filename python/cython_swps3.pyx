cimport cython_swps3

import numpy as np
cimport numpy as np
import json

#This function expects "normalized" byte arrays as input parameters, which means that 'A' should be normalized to 0
# 'B' to 1, 'C' to 2 etc... and by 0, 1, 2 NOT character '0', '1', '2' is meant
cpdef double swps3_alignScalar(np.ndarray[np.double_t,ndim=2] matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    return cython_swps3.python_alignScalar(<double*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


cpdef double swps3_alignByte(np.ndarray[np.int8_t,ndim=2] matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    return cython_swps3.python_alignByte(<signed char*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


cpdef double swps3_alignShort(np.ndarray[np.int8_t,ndim=2] matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold):
    return cython_swps3.python_alignShort(<signed char*> matrix.data, s1, ls1, s2, ls2, gapOpen, gapExt, threshold)


def alignScalar(s1, s2, env):
    return swps3_alignScalar(env.float64_Matrix, s1, len(s1), s2, len(s2), env.gapOpen, env.gapExt, env.threshold)


def normlaizedAlignScalar(s1, s2, env):
    return alignScalar(normalizeString(s1), normalizeString(s2), env)


def alignByte(s1, s2, env):
    return swps3_alignByte(env.int8_Matrix, s1, len(s1), s2, len(s2), env.gapOpen, env.gapExt, env.threshold)


def normalizedAlignByte(s1, s2, env):
    return alignByte(normalizeString(s1), normalizeString(s2), env)


def alignShort(s1, s2, env):
    return swps3_alignShort(env.int8_Matrix, s1, len(s1), s2, len(s2), env.gapOpen, env.gapExt, env.threshold)


def normalizedAlignShort(s1, s2, env):
    return alignShort(normalizeString(s1), normalizeString(s2), env)


def normalizeString(s):
    ret = ""
    for c in s:
        ret+=chr(ord(c) - ord('A'))

    return ret

class AlignmentEnvironment:
    def __init__(self):
        self.columns = "ARNDCQEGHILKMFPSTWYV"
        self.gapOpen = 0.0
        self.gapExt = 0.0
        self.threshold = 120.0
        self.float64_Matrix = np.ndarray(shape=(26, 26), dtype=np.float64)
        self.int8_Matrix = np.ndarray(shape=(26, 26), dtype=np.int8)

    def readFromFile(self, file):
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
        n = 26
        dtype = np.dtype(np.float64)
        nbytes = n * n * dtype.itemsize
        buf = np.empty(nbytes + 16, dtype=np.uint8)
        start_index = -buf.ctypes.data % 16
        a = buf[start_index:start_index + nbytes].view(dtype).reshape(n, n)

        #copy the not aligned values to the aligned NumPy matrix
        np.copyto(a, self.float64_Matrix)
        self.float64_Matrix = a
        self.int8_Matrix = self.float64_Matrix.astype(np.int8)