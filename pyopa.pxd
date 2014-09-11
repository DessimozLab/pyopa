cdef extern from "Python_extension.h":
    ctypedef struct ProfileByte:
        pass
    ctypedef struct ProfileShort:
        pass

    double c_align_scalar_reference_local(double* matrix, const char *s1, int ls1, const char *s2, int ls2,
                               double gap_open, double gap_ext, double threshold)
    double c_align_profile_byte_sse_local(ProfileByte* profile, const char *s2, int ls2,
                                         double gap_open, double gap_ext, double threshold)
    double c_align_profile_short_sse_local(ProfileShort* profile, const char *s2, int ls2,
                                          double gap_open, double gap_ext, double threshold)

    ProfileByte* c_create_profile_byte_sse_local(const char* query, int query_len, signed char* matrix)
    ProfileShort* c_create_profile_short_sse_local(const char* query, int query_len, signed short* matrix)

    void c_free_profile_byte_sse_local(ProfileByte* profile)
    void c_free_profile_short_sse_local(ProfileShort* profile)

    double c_align_double_local(double* matrix, const char *s1, int ls1, const char *s2, int ls2, double gap_open,
                          double gap_ext, double threshold, int* max1, int* max2)
    double c_align_double_global(double* matrix, const char *s1, int ls1, const char *s2, int ls2,
                                 double gap_open, double gap_ext)
    int c_align_strings(double* matrix, char *s1, int len1, char *s2, int len2, double escore, char *o1,
                        char *o2, double maxerr, double gap_open, double gap_ext)


cdef extern from "EstimatePam.h":
    ctypedef struct DayMatrix:
        pass

    DayMatrix* createDayMatrices(double* gapOpen, double* gapExt,
        double* pamDistances, long long* matrix_pointers, int DMSLen)
    void EstimatePam(char* o1, char* o2, int len, DayMatrix* DMS, int DMSLen,
            double* logPAM1, double* result)
    void freeDayMatrices(DayMatrix* DMS, int DMSLen)
    void CreateOrigDayMatrix(double* log_pam1, double pam, double* new_matrix)


cdef extern from "DynProgr_sse_double.h":
    ctypedef struct ProfileDouble:
        pass

    double align_double_local(ProfileDouble* profileDouble, const char *s2, int ls2, double gap_open,
        double gap_ext, double threshold, int* max1, int* max2)
    ProfileDouble* createProfileDoubleSSE(const char* s1, int ls1, double* matrix)
    void free_profile_double_sse(ProfileDouble* profile)