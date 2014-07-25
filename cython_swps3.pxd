cdef extern from "Python_extension.h":
    ctypedef struct ProfileByte:
        pass
    ctypedef struct ProfileShort:
        pass

    double c_align_scalar(double* matrix, const char *s1, int ls1, const char *s2, int ls2,
                               double gap_open, double ga_ext, double threshold)
    double c_align_profile_byte_sse(ProfileByte* profile, const char *s2, int ls2,
                                         double gap_open, double gap_ext, double threshold)
    double c_align_profile_short_sse(ProfileShort* profile, const char *s2, int ls2,
                                          double gap_open, double gap_ext, double threshold)

    ProfileByte* c_create_profile_byte_sse(const char* query, int query_len, signed char* matrix)
    ProfileShort* c_create_profile_short_sse(const char* query, int query_len, signed short* matrix)

    void c_free_profile_byte_sse(ProfileByte* profile)
    void c_free_profile_short_sse(ProfileShort* profile)