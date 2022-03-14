#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "htsutils.h"

// return
// -1: no seq found in file, not fasta or fastq
// 0: fasta
// 1: fastq
int is_fastq(char *fn) {
	int is_fq = -1;
	gzFile infp = xzopen(fn, "r");
    kseq_t *seq = kseq_init(infp);

    while (kseq_read(seq) >= 0) {
    	if (seq->qual.l > 0) is_fq = 1;
    	else is_fq = 0;
    	break;
    }
    err_gzclose(infp);
    kseq_destroy(seq);
    return is_fq;
}
