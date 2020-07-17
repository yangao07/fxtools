#ifndef __PARSE_H__
#define __PARSE_H__
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <ctype.h>
#include "utils.h"
#include "kstring.h"
#include "htslib/sam.h"

static inline int aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

static long int _strtol10(const char *str, char **endptr) {
    long int res = 0; unsigned d;
    char *s;
    for (s = (char*)str, d = s[0]-'0'; d < 10; ++s, d = s[0]-'0')
        res = res * 10 + d;
    if (endptr) *endptr = s;
    return res;
}

static unsigned long int _strtoul10(const char *str, char **endptr) {
    unsigned long int res = 0; unsigned d;
    char *s;
    for (s = (char*)str, d = s[0]-'0'; d < 10; ++s, d = s[0]-'0')
        res = res * 10 + d;
    if (endptr) *endptr = s;
    return res;
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += aux_type2size(type); \
	} while(0)

uint8_t *aux_get(int l_data, const uint8_t *data, const char tag[2])
{
	const uint8_t *s = data;
	int y = tag[0]<<8 | tag[1];
	while (s < data + l_data) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return (uint8_t*)s;
		__skip_tag(s);
	}
	return 0;
}

// s MUST BE returned by aux_get()
int aux_del(int l_data, uint8_t *data, uint8_t *s)
{
	uint8_t *p;
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, l_data - (s - data));
	return l_data - (s - p);
}

int aux_parse(char *s, uint8_t **data, int *max)
{
	char *q, *p;
	kstring_t str;
	if (s == 0) return 0;
	str.l = 0, str.m = *max, str.s = (char*)*data;
	if (*s == '\t') ++s;
	for (p = q = s;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'I' || q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
				int type = q[3];
				kputsn_(q, 2, &str);
				q += 5;
				if (type == 'A') {
					kputc_('A', &str);
					kputc_(*q, &str);
                } else if (type == 'I') {
					uint32_t x;
					// x = strtol(q, &q, 10);
                    x = _strtoul10(q, &q);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'i') {
					int32_t x;
					// x = strtol(q, &q, 10);
                    x = _strtol10(q, &q);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'f') {
					float x;
					x = strtod(q, &q);
					kputc_('f', &str); kputsn_(&x, 4, &str);
				} else if (type == 'Z') {
					kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
				} else if (type == 'B') {
					type = *q++; // q points to the first ',' following the typing byte
					if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
						int32_t n;
						char *r;
						for (r = q, n = 0; *r; ++r)
							if (*r == ',') ++n;
						kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
						// TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
						if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						// else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'i') while (q + 1 < p) { int32_t  x = _strtol10(q + 1, &q); kputsn_(&x, 4, &str); }
						// else if (type == 'I') while (q + 1 < p) { uint32_t x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'I') while (q + 1 < p) { uint32_t x = _strtoul10(q + 1, &q); kputsn_(&x, 4, &str); }
						else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
					}
				} // should not be here, as we have tested all types
			}
			q = p + 1;
			if (c == 0) break;
		}
	}
	if (str.l > 0 && str.l == str.m) ks_resize(&str, str.l + 1);
	if (str.s) str.s[str.l] = 0;
	*max = str.m, *data = (uint8_t*)str.s;
	return str.l;
}

int _bam_cigar2qlen(int n_cigar, uint32_t *cigar) {
    int i, len = 0, l, op;
    for (i = 0; i < n_cigar; ++i) {
        l = cigar[i] >> 4; 
        op = cigar[i] & 0xf;
        if (op != BAM_CDEL && op != BAM_CREF_SKIP)
            len += l;
    }
    return len;
}

int cigar_str_parse(char *cigar_str, int *match, int *mis, int *ins, int *del, int *skip, int *clip) {
    int equal = 0, diff = 0;
    *match = *mis = *ins = *del = *skip = *clip = 0;
    int len;
    char *p= cigar_str, *end;
    for (; *p; ) {
        len = strtol(p, &end, 10);
        switch (*end) {
            case 'M': *match += len; break;
            case '=': equal += len; break;
            case 'X': diff += len; break;
            case 'I': *ins += len; break;
            case 'D': *del += len; break;
            case 'N': *skip += len; break;
            case 'S': *clip += len; break;
            case 'H': *clip += len; break;
            default: err_fatal(__func__, "Cigar string Error. (%s)\n", cigar_str);
        }
        p = end+1;
    }

    if (equal != 0 || diff != 0) {
        *mis = diff; *match = equal;
    }
    return 0;
}

int cigar_parse(bam1_t *b, uint32_t *cigar, int cigar_len, int *match, int *mis, int *ins, int *del, int *skip, int *clip) {
    int i, md, equal = 0, diff = 0;
    *match = *mis = *ins = *del = *skip = *clip = 0;
    for (i = 0; i < cigar_len; ++i) {
        uint32_t c = cigar[i];
        int len = bam_cigar_oplen(c);
        switch (bam_cigar_op(c)) {
            case BAM_CMATCH: *match += len; break;
            case BAM_CEQUAL: equal += len; break;
            case BAM_CDIFF: diff += len; break;
            case BAM_CINS: *ins += len; break;
            case BAM_CDEL: *del += len; break;
            case BAM_CREF_SKIP: *skip += len; break;
            case BAM_CSOFT_CLIP: *clip += len; break;
            case BAM_CHARD_CLIP: *clip += len; break;
            default : err_fatal_simple("Cigar ERROR.\n");
        }
    }
    if (equal == 0 && diff == 0) {
        uint8_t *p = bam_aux_get(b, "NM");
        if (p == 0) p = bam_aux_get(b, "nM");
        if (p == 0) {
            err_fatal_core(__func__, "%s No \"NM\" tag.\n", bam_get_qname(b));
            return 0;
        }
        md = bam_aux2i(p);
        *mis = md - *ins - *del;
        *match = *match - *mis;
    } else {
        *mis = diff; *match = equal;
    }
    return 0;
}

#endif
