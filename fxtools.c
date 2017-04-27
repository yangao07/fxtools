#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdint.h>
#include "utils.h"
#include "fxtools.h"
#include "kseq.h"
#include "htslib/htslib/sam.h"

#define _ll_t long long

KSEQ_INIT(gzFile, gzread)

int usage(void)
{
    fprintf(stderr, "Program: fxtools (fasta and fastq data tools)\n");
    fprintf(stderr, "Usage:   fxtools <command> [options]\n\n");
    fprintf(stderr, "Command: filter (fl)         filter fa/fq sequences with specified length bound.\n");
    fprintf(stderr, "         fq2fa (qa)          convert FASTQ format data to FASTA format data.\n");
    fprintf(stderr, "         fa2fq (aq)          convert FASTA format data to FASTQ format data.\n");
    fprintf(stderr, "         re-co (rc)          convert DNA sequence(fa/fq) to its reverse-complementary sequence.\n");
    fprintf(stderr, "         seq-display (sd)    display a specified region of FASTA/FASTQ file.\n");
    fprintf(stderr, "         cigar-parse (cp)    parse the given cigar(stdout).\n");
    fprintf(stderr, "         length-parse (lp)   parse the length of sequences in fa/fq file.\n");
    fprintf(stderr, "         merge-fa (mf)       merge the reads with same read name in fasta/fastq file.\n");
    fprintf(stderr, "         error-parse (ep)    parse indel and mismatch error based on CIGAR and NM in bam file.\n");
    //fprintf(stderr, "      ./fa_filter in.fa out.fa low-bound upper-bound(-1 for no bound)\n");
    fprintf(stderr, "\n");
    return 1;
}

int fxt_filter(int argc, char* argv[])
{
    if (argc != 4) 
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools filter <in.fa/fq> <lower-bound> <upper-bound>(-1 for NO bound)\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp;
    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) infp = gzdopen(fileno(stdin), "r");
    else infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_filter] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq;
    seq = kseq_init(infp);
    int64_t low = atoi(argv[2]);
    int64_t upper = atoi(argv[3]);
    while (kseq_read(seq) >= 0)
    {
        if ((low != -1 && (int64_t)seq->seq.l < low) || (upper != -1 && (int64_t)seq->seq.l > upper))
            continue;
        if (seq->qual.l != 0)
        {
            fprintf(stdout, "@%s\n", seq->name.s);
            fprintf(stdout, "%s\n", seq->seq.s);
            fprintf(stdout, "+\n");
            fprintf(stdout, "%s\n", seq->qual.s);
        }
        else
        {
            fprintf(stdout, ">%s\n", seq->name.s);
            fprintf(stdout, "%s\n", seq->seq.s);
        }
    }

    gzclose(infp);
    return 0;
}

int fxt_fq2fa(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools fq2fa <in.fq> <out.fa>\n\n");
        exit(-1);
    } 
    gzFile infp;
    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) infp = gzdopen(fileno(stdin), "r");
    else infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_fq2fa] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = fopen(argv[2], "w");

    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, ">%s\n", seq->name.s);
        fprintf(outfp, "%s\n", seq->seq.s);
    }

    gzclose(infp);
    fclose(outfp);
    return 0;
}

int fxt_fa2fq(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools fa2fq <in.fa> <out.fq>\n\n");
        exit(-1);
    } 
    gzFile infp;
    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) infp = gzdopen(fileno(stdin), "r");
    else infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_fa2fq] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = fopen(argv[2], "w");

    int64_t i;
    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, "@%s\n", seq->name.s);
        fprintf(outfp, "%s\n", seq->seq.s);
        fprintf(outfp, "+\n");
        for (i = 0; i < (int64_t)seq->seq.l; ++i) fprintf(outfp, "!");
        fprintf(outfp, "\n");
    }

    gzclose(infp);
    fclose(outfp);
    return 0;
}

int fxt_re_co(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools re-co in.fa/fq out.fa\n"); fprintf(stderr, "\n");
        return 1;
    }
    gzFile readfp;
    kseq_t *read_seq;
    int seq_len = 100000;
    char *seq = (char*)malloc(seq_len*sizeof(char));
    int8_t *seq_n = (int8_t*)malloc(seq_len*sizeof(int8_t));
    FILE *out = fopen(argv[2], "w");

    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) readfp = gzdopen(fileno(stdin), "r");
    else readfp = gzopen(argv[1], "r");
    read_seq = kseq_init(readfp);

    int len, i;

    while ((len = kseq_read(read_seq)) > 0)
    {
       if (len > seq_len)
       {
           seq_len <<= 1;
           seq = (char*)realloc(seq, seq_len*sizeof(char));
           seq_n = (int8_t*)realloc(seq_n, seq_len*sizeof(char));
           if (seq == NULL || seq_n == NULL)
           {
               fprintf(stderr, "memory is not enough.\n");
               exit(-1);
           }
       } 
       seq = read_seq->seq.s;
       for (i = 0; i < len; i++)
           seq_n[i] = nt_table[(int)seq[i]];
       fprintf(out, ">%s_re-co:\n", read_seq->name.s);
       for (i = len - 1; i>=0; i--)
       {
           if (seq_n[i] != 4) 
               fprintf(out, "%c", nt_char[3-(int)seq_n[i]]);
            else
                fprintf(out, "N");
       }
       fprintf(out, "\n");
    }

    gzclose(readfp);
    kseq_destroy(read_seq);
    fclose(out);

    return 0;
}

int fxt_seq_dis(int argc, char *argv[])
{
    if (argc != 5)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools seq-display <in.fa/fq> chr/read_name start_pos(1-based) length\n");
        fprintf(stderr, "\n"); 
        exit(-1);
    }
    gzFile infp;
    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) infp = gzdopen(fileno(stdin), "r");
    else infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_seq_dis] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    char name[1024]; uint64_t start; int len;
    strcpy(name, argv[2]); start = atol(argv[3]); len = atoi(argv[4]);
    kseq_t *seq;
    seq = kseq_init(infp);
    while (kseq_read(seq) >= 0)
    {
        if (strcmp(seq->name.s, name) == 0) {
            if (start > seq->seq.l) {
                fprintf(stderr, "[fxt_seq_dis] Error: START_POS is longger than the length of chr/read. (%lld > %lld)\n", (_ll_t)start, (_ll_t)seq->seq.l);
                exit(-1);
            } else if (start + len - 1 > seq->seq.l) {
                fprintf(stderr, "[fxt_seq_dis] Error: START_POS+LEN is longger than the length of chr/read. (%lld > %lld)\n", (_ll_t)start+len-1, (_ll_t)seq->seq.l);
            } else {
                seq->seq.s[start+len-1] = '\0';
            }
            fprintf(stdout, "%s\n", seq->seq.s+start-1);
            gzclose(infp);
            return 0;
        }
    }
    gzclose(infp);
    fprintf(stderr, "[fxt_seq_dis] Error: Can't find %s in %s.\n", name, argv[1]);
    return 0;
}

int fxt_cigar_parse(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools cigar-parse <input-cigar>\n\n");
        return 1;
    }
    int cigar_len, i, seq_len, ref_len;
    int c;
    long x, op[11] = {0};
    char *s, *t;

    cigar_len = seq_len = ref_len = 0;
    for (s = argv[1]; *s; )
    {
        x = strtol(s, &t, 10);  
        /*if (x == 0)
        {
            fprintf(stderr, "%s\n",s);
            fprintf(stderr, "[fxtools cigar-parse] Cigar ERROR 1.\n");
            exit(-1);
        }*/
        c = toupper(*t);
        switch (c)
        {
            case 'M':   op[CMATCH]+=x, seq_len+=x, ref_len+=x;    break;
            case 'I':   op[CINS]+=x, seq_len+=x;      break;
            case 'D':   op[CDEL]+=x, ref_len+=x;      break;
            case 'N':   op[CREF_SKIP]+=x, ref_len+=x;     break;
            case 'S':   op[CSOFT_CLIP]+=x, seq_len+=x;    break;
            case 'H':   op[CHARD_CLIP]+=x;    break;
            case 'P':   op[CPAD]+=x;          break;
            case '=':   op[CEQUAL]+=x, seq_len+=x, ref_len+=x;    break;
            case 'X':   op[CDIFF]+=x, seq_len+=x, ref_len+=x; break;
            case 'B':   op[CBACK]+=x; break;  
			case 'V':	op[CINV]+=x, seq_len+=x, ref_len+=x;	break;
            default:    fprintf(stderr, "[fxtools cigar-parse] Cigar ERROR 2.\n"); exit(-1); break;
        }
        //modify variable directly OR use a auxiliary-variable
        ++cigar_len;
        s = t+1;
    }
    for (i = 0; i < 11; ++i)
    {
        if (op[i] != 0) fprintf(stdout, "%ld%c\t", op[i], CIGAR_STR[i]);
    }
    fprintf(stdout, "\nseq-len: %d\nref-len: %d\n", seq_len, ref_len);
    return 0;
}

int fxt_len_parse(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools length-parse <in.fa/fq>\n");
        fprintf(stderr, "\n"); 
        exit(-1);
    }
    gzFile infp;
    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) infp = gzdopen(fileno(stdin), "r");
    else infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_len_parse] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq;
    seq = kseq_init(infp);
    while (kseq_read(seq) >= 0)
    {
        fprintf(stdout, "%s\t%d\n", seq->name.s, (int)seq->seq.l);       
    }

    gzclose(infp);
    return 0;
}

int fxt_merge_fa(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools merge_fa <in.fa/fq> <out.fa/fq>\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp;
    if (strcmp(argv[1],"-") == 0 || strcmp(argv[1], "stdin") == 0) infp = gzdopen(fileno(stdin), "r");
    else infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_merge_fa] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq = kseq_init(infp);
    FILE *outfp = fopen(argv[2], "w");
    char read_name[1024]; 
    char *read_seq = (char*)calloc(10, 1);
    char *read_qual = (char*)calloc(10, 1);
    int w_seq_n=0;

    while (kseq_read(seq) >= 0)
    {
        if (strcmp(seq->name.s, read_name) == 0) {
            w_seq_n += seq->seq.l;
            read_seq = (char*)realloc(read_seq, w_seq_n+seq->seq.l);
            strcat(read_seq, seq->seq.s);
            read_seq[w_seq_n] = 0;
            if (seq->qual.l > 0) {
                read_qual = (char*)realloc(read_qual, w_seq_n+seq->seq.l);
                strcat(read_qual, seq->qual.s);
                read_qual[w_seq_n] = 0;
            }
        } else {
            if (w_seq_n > 0) {
                if (seq->qual.l > 0) { // fastq
                    fprintf(outfp, "@%s\n", read_name);
                    fprintf(outfp, "%s\n", read_seq);
                    fprintf(outfp, "+\n");
                    fprintf(outfp, "%s\n", read_qual);
                } else { // fasta
                    fprintf(outfp, ">%s\n", read_name);
                    fprintf(outfp, "%s\n", read_seq);
                }
            }
            w_seq_n = seq->seq.l;
            strcpy(read_name, seq->name.s);
            read_seq = (char*)realloc(read_seq, w_seq_n+seq->seq.l);
            strcpy(read_seq, seq->seq.s);
            read_seq[w_seq_n] = 0;
            if (seq->qual.l > 0) {
                read_qual = (char*)realloc(read_qual, w_seq_n+seq->seq.l);
                strcpy(read_qual, seq->qual.s);
                read_qual[w_seq_n] = 0;
            }
        }
    }
    if (w_seq_n > 0) { // last read
        if (seq->qual.l > 0) { // fastq
            fprintf(outfp, "@%s\n", read_name);
            fprintf(outfp, "%s\n", read_seq);
            fprintf(outfp, "+\n");
            fprintf(outfp, "%s\n", read_qual);
        } else { // fasta
            fprintf(outfp, ">%s\n", read_name);
            fprintf(outfp, "%s\n", read_seq);
        }
    }
    gzclose(infp);
    fclose(outfp);
    return 0;
}

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

int fxt_error_parse(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools error-parse <input.bam> > error.out\n\n");
        return 1;
    }
    fprintf(stdout, "READ_NAME\tREAD_LEN\tINS\tDEL\tMIS\tMATCH\tCLIP\tSKIP\n");
    long long tol_len=0, tol_ins=0, tol_del=0, tol_mis=0, tol_match=0, tol_clip=0, tol_skip=0;
    int i, seq_len, md, ins, del, mis, match, clip, skip;

    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1(); 

    while (sam_read1(in, h, b) >= 0) {
        seq_len = b->core.l_qseq;
        md = 0, ins = 0, del = 0, mis = 0, match = 0, clip = 0, skip = 0;
        if (!bam_unmap(b)) {
            uint32_t *cigar = bam_get_cigar(b); int cigar_len = b->core.n_cigar;
            for (i = 0; i < cigar_len; ++i) {
                uint32_t c = cigar[i];
                int len = bam_cigar_oplen(c);
                switch (bam_cigar_op(c)) {
                    case BAM_CMATCH: match += len; break;
                    case BAM_CINS: ins += len; break;
                    case BAM_CDEL: del += len; break;
                    case BAM_CREF_SKIP: skip += len; break;
                    case BAM_CSOFT_CLIP: clip += len; break;
                    case BAM_CHARD_CLIP: clip += len; break;
                    default : err_fatal_simple("Cigar ERROR.\n");
                }
            }
            uint8_t *p = bam_aux_get(b, "NM");
            if (p == 0) p = bam_aux_get(b, "nM");
            if (p == 0) {
                err_fatal_core(__func__, "%s No \"NM\" tag.\n", bam_get_qname(b));
                return 0;
            }
            md = bam_aux2i(p);
            mis = md - ins - del;
            match = match - mis;
        }
        tol_len += seq_len; tol_ins += ins; tol_del += del; tol_mis += mis; tol_match += match; tol_clip += clip; tol_skip += skip;

        fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", bam_get_qname(b), seq_len, ins, del, mis, match, clip, skip);
    }
    fprintf(stdout, "%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n", "Total", tol_len, tol_ins, tol_del, tol_mis, tol_match, tol_clip, tol_skip);
    return 0;
}

int main(int argc, char*argv[])
{
    if (argc < 2) return usage();
    if (strcmp(argv[1], "filter") == 0 || strcmp(argv[1], "fl") == 0) fxt_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "fq2fa") == 0 || strcmp(argv[1], "qa") == 0) fxt_fq2fa(argc-1, argv+1);
    else if (strcmp(argv[1], "fa2fq") == 0 || strcmp(argv[1], "aq") == 0) fxt_fa2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "re-co") == 0 || strcmp(argv[1], "rc") == 0) fxt_re_co(argc-1, argv+1);
    else if (strcmp(argv[1], "seq-display") == 0 || strcmp(argv[1], "sd") == 0) fxt_seq_dis(argc-1, argv+1);
    else if (strcmp(argv[1], "cigar-parse") == 0 || strcmp(argv[1], "cp") == 0) fxt_cigar_parse(argc-1, argv+1);
    else if (strcmp(argv[1], "length-parse") == 0 || strcmp(argv[1], "lp") == 0) fxt_len_parse(argc-1, argv+1);
    else if (strcmp(argv[1], "merge-fa") == 0 || strcmp(argv[1], "mf") == 0) fxt_merge_fa(argc-1, argv+1);
    else if (strcmp(argv[1], "error-parse") == 0 || strcmp(argv[1], "ep") == 0) fxt_error_parse(argc-1, argv+1);
    else {fprintf(stderr, "unknow command [%s].\n", argv[1]); return 1; }

    return 0;
}
