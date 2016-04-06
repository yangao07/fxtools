#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "fxtools.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

int usage(void)
{
    fprintf(stderr, "Program: fxtools (fasta and fastq data tools)\n");
    fprintf(stderr, "Usage:   fxtools <command> [options]\n\n");
    fprintf(stderr, "Command: filter            filter fa/fq sequences with specified length bound.\n");
    fprintf(stderr, "         fq2fa             convert FASTQ format data to FASTA format data.\n");
    fprintf(stderr, "         fa2fq             convert FASTA format data to FASTQ format data.\n");
    fprintf(stderr, "         re-co             convert DNA sequence(fa/fq) to its reverse-complementary sequence.\n");
    fprintf(stderr, "         cigar-parse       parse the given cigar(stdout).\n");
    fprintf(stderr, "         length-parse      parse the length of sequences in fa/fq file.\n");
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
    gzFile infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_filter] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq;
    seq = kseq_init(infp);
    int low = atoi(argv[2]);
    int upper = atoi(argv[3]);
    while (kseq_read(seq) >= 0)
    {
        if ((low != -1 && seq->seq.l < low) || (upper != -1 && seq->seq.l > upper))
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
    gzFile infp = gzopen(argv[1], "r");
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
    gzFile infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_fa2fq] Can't open %s.\n", argv[1]);
        exit(-1);
    }
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = fopen(argv[2], "w");

    int i;
    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, "@%s\n", seq->name.s);
        fprintf(outfp, "%s\n", seq->seq.s);
        fprintf(outfp, "+\n");
        for (i = 0; i < seq->seq.l; ++i) fprintf(outfp, "!");
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

    readfp = gzopen(argv[1], "r");
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

int fxt_cigar_parse(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools cigar-parse <input-cigar>\n\n");
        return 1;
    }
    int cigar_len, i, seq_len;
    int c;
    long x, op[11] = {0};
    char *s, *t;

    cigar_len = seq_len = 0;
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
            case 'M':   op[CMATCH]+=x, seq_len+=x;    break;
            case 'I':   op[CINS]+=x, seq_len+=x;      break;
            case 'D':   op[CDEL]+=x;      break;
            case 'N':   op[CREF_SKIP]+=x;     break;
            case 'S':   op[CSOFT_CLIP]+=x, seq_len+=x;    break;
            case 'H':   op[CHARD_CLIP]+=x;    break;
            case 'P':   op[CPAD]+=x;          break;
            case '=':   op[CEQUAL]+=x, seq_len+=x;    break;
            case 'X':   op[CDIFF]+=x, seq_len+=x; break;
            case 'B':   op[CBACK]+=x; break;  
			case 'V':	op[CINV]+=x, seq_len+=x;	break;
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
    fprintf(stdout, "\nseq-len: %d\n", seq_len);
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
    gzFile infp = gzopen(argv[1], "r");
    if (infp == NULL)
    {
        fprintf(stderr, "[fxt_filter] Can't open %s.\n", argv[1]);
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

int main(int argc, char*argv[])
{
    if (argc < 2) return usage();
    if (strcmp(argv[1], "filter") == 0) fxt_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "fq2fa") == 0) fxt_fq2fa(argc-1, argv+1);
    else if (strcmp(argv[1], "fa2fq") == 0) fxt_fa2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "re-co") == 0) fxt_re_co(argc-1, argv+1);
    else if (strcmp(argv[1], "cigar-parse") == 0) fxt_cigar_parse(argc-1, argv+1);
    else if (strcmp(argv[1], "length-parse") == 0) fxt_len_parse(argc-1, argv+1);
    else {fprintf(stderr, "unknow command [%s].\n", argv[1]); return 1; }

    return 0;
}
