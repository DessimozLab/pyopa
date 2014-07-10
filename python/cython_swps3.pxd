cdef extern from "Python_extension.h":
    double python_alignScalar(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold)
    double python_alignByte(signed char* matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold)
    double python_alignShort(signed char* matrix, const char *s1, int ls1, const char *s2, int ls2, double gapOpen, double gapExt, double threshold)