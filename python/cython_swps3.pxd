cdef extern from "Python_extension.h":
    ctypedef struct ProfileByte:
        pass
    ctypedef struct ProfileShort:
        pass

    double python_alignScalar(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold)
    double python_alignByteProfileSSE(ProfileByte* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold)
    double python_alignShortProfileSSE(ProfileShort* profile, const char *s2, int ls2, double gapOpen, double gapExt, double threshold)

    ProfileByte* python_createByteProfileSSE(const char* query, int queryLen, signed char* matrix)
    ProfileShort* python_createShortProfileSSE(const char* query, int queryLen, signed short* matrix)

    void python_freeProfileByteSSE(ProfileByte* profile)
    void python_freeProfileShortSSE(ProfileShort* profile)