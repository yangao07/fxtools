#ifndef HTSUTILS_H
#define HTSUTILS_H

#include <zlib.h>
#include "utils.h"
#include "kseq.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"

KSEQ_INIT(gzFile, gzread)

#ifdef __cplusplus
extern "C" {
#endif

int is_fastq(char *fn);

#ifdef __cplusplus
}
#endif

#endif
