#ifndef __SBCD_WL_H
#define __SBCD_WL_H

#include "seqmatch.h"
#include "ErrorWarning.h"

extern "C" {
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/kseq.h"
#include "htslib/htslib/hts.h"
#include "htslib/htslib/tbx.h"
}

class SbcdWL {

private:
    Parameters& P; //reference to the parameters, will be initialized on construction
    htsFile* hp;
    kstring_t str;         // kstring_t object to store the string
    int32_t* fields;       // indices of the starting point to each field
    int32_t lstr, nfields, delimiter;

public:

    BcdMatch<uint64_t> bcd_match;
    std::unordered_map<uint64_t, uint64_t> val_idx_map;
    uint64_t n_uniq_val;
	int32_t bcd_len, kmer_size;
    int32_t max_mismatch = 1;
	std::string crd_tag, bcd_tag, bcd_tag_org;
    bool active = false;
    bool unique_value_map = false;
    char bcd_tag_c[2], crd_tag_c[2], bcd_tag_org_c[2];

    SbcdWL(Parameters& _P) : P(_P) {}
    SbcdWL(Parameters& _P, int32_t _bcL, int32_t _k, bool _exact = false, bool _ar = true, bool _aq = true) : P(_P) {
        init(_bcL, _k, _exact, _ar, _aq);
    }
    void init(int32_t _bcL, int32_t _k, bool _exact = false, bool _ar = true, bool _aq = true) {
        bcd_len = _bcL;
        kmer_size = _k;
        bcd_match.init(bcd_len, kmer_size, _exact, _ar, _aq);
        str.l = str.m = 0; str.s = NULL;
        delimiter = 0;
        nfields = 0;
        fields = NULL;
        active = true;
    }
    int32_t loadWL(const char* filename, uint32_t _icol_s, uint32_t _icol_x, uint32_t _icol_y, bool _list_val = false);
    int32_t query(const char* seq, char* cb, int32_t& x, int32_t& y);

};

#endif // __SBCD_WL_H
