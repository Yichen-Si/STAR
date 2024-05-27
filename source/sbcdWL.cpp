#include "sbcdWL.h"

int32_t SbcdWL::loadWL(const char* filename, uint32_t icol_s, uint32_t icol_x, uint32_t icol_y, bool _list_val) {
    if (!active) {
        return 0;
    }
    unique_value_map = _list_val;
    if (icol_s == icol_x || icol_s == icol_y || icol_x == icol_y) {
        exitWithError("Error: Column indices must be unique", std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    hp = hts_open(filename, "r");
    if (hp == NULL) {
        exitWithError("Error: Cannot open the white list (" + std::string(filename) + ")", std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    }
    int32_t min_fields = std::max(icol_s, std::max(icol_x, icol_y)) + 1;
    lstr = hts_getline(hp, KS_SEP_LINE, &str);
    fields = ksplit(&str, delimiter, &nfields);
    if (nfields < min_fields) {
        exitWithError("Error: Column indices must be from 0 to ncol - 1", std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    uint64_t ref_tot = 0, ref_ambig = 0, ref_skip = 0;
    n_uniq_val = 0;
    while(true) {
        ref_tot++;
        if (ref_tot % 500000 == 0) {
            printf("SbcdWL::loadWL: processed %lu records, %lu contain one ambiguous base, %lu are skipped\n", ref_tot, ref_ambig, ref_skip);
        }
        uint64_t xy = (uint64_t) atoi(&str.s[fields[icol_x]]) << 32 | (uint64_t) atoi(&str.s[fields[icol_y]]);
        int32_t ret = bcd_match.add_ref(&str.s[fields[icol_s]], xy);
        if (ret < 0) {
            ref_skip++;
        } else if (ret > 0) {
            ref_ambig++;
        }
        if (unique_value_map && ret >= 0) {
            auto it = val_idx_map.emplace(xy, n_uniq_val);
            if (it.second) {
                n_uniq_val++;
            }
        }
        lstr = hts_getline(hp, KS_SEP_LINE, &str);
        if (lstr <= 0) {
            break;
        }
        if ( fields != NULL ) { free(fields); fields = NULL; }
        fields = ksplit(&str, delimiter, &nfields);
        if (nfields < min_fields) {
            exitWithError("Error: ill-formated line in the white list", std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        if (P.debug % 3 == 1 && ref_tot > 500000) {break;}
    }
    bcd_match.process_ambig_ref();
    P.inOut->logProgress << "Read white list with " << ref_tot << " record including " << ref_ambig << " with one ambiguous base, skipped " << ref_skip << std::endl << std::flush;
    return 1;
}


int32_t SbcdWL::query(const char* seq, char* cb, int32_t& x, int32_t& y) {
    if (!active) {
        exitWithError("Error: The white list is not initialized", std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
    }
    uint64_t nt4cb, xy;
    int32_t nmiss = bcd_match.query(seq, nt4cb, xy);
    if (nmiss < 0) {
        return nmiss;
    }
    x = (int32_t) (xy >> 32);
    y = (int32_t) (xy & 0xFFFFFFFF);
    if (nmiss > 0) {
        cb[bcd_len] = '\0';
        for (int32_t i = 0; i < bcd_len; i++) {
            cb[i] = seq_nt4_rev_table[nt4cb&0x3];
            nt4cb >>= 2;
        }
    }
    return nmiss;
}
