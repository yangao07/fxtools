
//copy from samtools-0.1.19
#define CIGAR_STR "MIDNSHP=XB"
/*
  CIGAR operations.
 */
/*! @abstract CIGAR: M = match or mismatch*/
#define CMATCH      0
/*! @abstract CIGAR: I = insertion to the reference */
#define CINS        1
/*! @abstract CIGAR: D = deletion from the reference */
#define CDEL        2
/*! @abstract CIGAR: N = skip on the reference (e.g. spliced alignment) */
#define CREF_SKIP   3
/*! @abstract CIGAR: S = clip on the read with clipped sequence
  present in qseq */
#define CSOFT_CLIP  4
/*! @abstract CIGAR: H = clip on the read with clipped sequence trimmed off */
#define CHARD_CLIP  5
/*! @abstract CIGAR: P = padding */
#define CPAD        6
/*! @abstract CIGAR: equals = match */
#define CEQUAL      7
/*! @abstract CIGAR: X = mismatch */
#define CDIFF       8
#define CBACK       9

static int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
	};

char nt_char[5] = { 'A', 'C', 'G', 'T', 'N' };
