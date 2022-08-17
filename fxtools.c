#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <stdint.h>
#include <locale.h>
#include "fxtools.h"
#include "parse.h"
#include "khash.h"
#include "htsutils.h"

#define _ll_t long long

KHASH_MAP_INIT_STR(str, uint32_t)

char bam_qseq_table[17] = "NACNGNNNTNNNNNNN";
char bam_rc_qseq_table[17] = "NTGNCNNNANNNNNNN";

unsigned char nt16_bam_qseq_table[16] = {
    4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4
};
unsigned char nt16_bam_qseq_rc_table[16] = {
    4, 3, 2, 4, 1, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4
};

int usage(void)
{
    fprintf(stderr, "Program: fxtools (light-weight processing tool for FASTA, FASTQ and BAM format data)\n");
    fprintf(stderr, "Usage:   fxtools <command> [options]\n\n");
    fprintf(stderr, "Command: \n");
    fprintf(stderr, "         filter (fl)           filter fa/fq sequences with specified length boundary.\n");
    fprintf(stderr, "         filter-name (fn)      filter fa/fq sequences with specified name.\n");
    fprintf(stderr, "         filter-qual (fq)      filter fq sequences with specified min. average quality.\n");
    fprintf(stderr, "         filter-bam (fb)       filter bam/sam records with specified read length boundary.\n");
    fprintf(stderr, "         filter-bam-name (fbn) filter bam/sam records with specified read name.\n");
    fprintf(stderr, "         split-fx (sx)         split fa/fq file into multipule files.\n");
    fprintf(stderr, "         fq2fa (qa)            convert FASTQ format data to FASTA format data.\n");
    fprintf(stderr, "         fa2fq (aq)            convert FASTA format data to FASTQ format data.\n");
    fprintf(stderr, "         bam2bed (bb)          convert BAM file to BED file. seperated exon regions for spliced BAM\n");
    fprintf(stderr, "         bam2fx (bf)           convert BAM file to FASTA/FASTQ file.\n");
    fprintf(stderr, "         re-co (rc)            convert DNA sequence(fa/fq) to its reverse-complementary sequence.\n");
    fprintf(stderr, "         seq-display (sd)      display a specified region of FASTA/FASTQ file.\n");
    fprintf(stderr, "         cigar-parse (cp)      parse the given cigar(stdout).\n");
    fprintf(stderr, "         length-parse (lp)     parse the length of sequences in fa/fq file.\n");
    fprintf(stderr, "         merge-fa (mf)         merge the reads with same read name in fasta/fastq file.\n");
    fprintf(stderr, "         merge-filter-fa (mff) merge and filter the reads with same read name in fasta file.\n");
    fprintf(stderr, "         duplicate-seq (ds)    duplicate all the read sequences with specific copy number.\n");
    fprintf(stderr, "         duplicate-read (dd)   duplicate all the read records with specific copy number.\n");
    fprintf(stderr, "         error-parse (ep)      parse indel and mismatch error based on CIGAR and NM in SAM/BAM/GAF file.\n");
    fprintf(stderr, "         dna2rna (dr)          convert DNA fa/fq to RNA fa/fq.\n");
    fprintf(stderr, "         rna2dna (rd)          convert RNA fa/fq to DNA fa/fq.\n");
    fprintf(stderr, "         trim (tr)             trim polyA tail(polyT head).\n");
    fprintf(stderr, "         trimF (tf)            trim and filter with polyA tail(polyT head). Only polyA reads will be kept.\n");
    fprintf(stderr, "         bam-add-seq (bt)      add polyA sequence to BAM record.\n");
    fprintf(stderr, "\n");
    return 1;
}

void print_seq(FILE *out, kseq_t *seq) {
    if (seq->seq.l == 0) return;
    if (seq->qual.l != 0)
    {
        fprintf(out, "@%s", seq->name.s);
        if (seq->comment.l > 0) fprintf(out, " %s", seq->comment.s);
        fprintf(out, "\n");
        fprintf(out, "%s\n", seq->seq.s);
        fprintf(out, "+\n");
        fprintf(out, "%s\n", seq->qual.s);
    }
    else
    {
        fprintf(out, ">%s", seq->name.s);
        if (seq->comment.l > 0) fprintf(out, " %s", seq->comment.s);
        fprintf(out, "\n");
        fprintf(out, "%s\n", seq->seq.s);
    }
}

// return 1 if filename has suffix of _suf_
// else 0
int check_suf(char *filename, char suf[]) {
    int i;
    int nl = strlen(filename), sl = strlen(suf);
    if (nl < sl) return 0;
    for (i = 0; i < sl; ++i) {
        if (toupper(filename[nl-i-1]) != toupper(suf[sl-i-1])) return 0;
    }
    return 1;
}

int fxt_filter(int argc, char* argv[]) {
    if (argc != 4) 
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools filter <in.fa/fq> <lower-bound> <upper-bound>(-1 for NO bound) > <out.fa/fq>\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[1], "r");
    FILE *out = stdout;
    kseq_t *seq;
    seq = kseq_init(infp);
    int64_t low = atoi(argv[2]);
    int64_t upper = atoi(argv[3]);
    while (kseq_read(seq) >= 0)
    {
        if ((low != -1 && (int64_t)seq->seq.l < low) || (upper != -1 && (int64_t)seq->seq.l > upper))
            continue;
        print_seq(out, seq);
    }

    err_fclose(out);
    kseq_destroy(seq);
    err_gzclose(infp);
    return 0;
}

int fxt_filter_bam(int argc, char *argv[]) {
    if (argc != 4) 
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools filter-bam <in.bam/sam> <lower-bound> <upper-bound>(-1 for NO bound) > <out.bam>\n");
        fprintf(stderr, "\n");
        exit(-1);
    }

    samFile *in, *out; bam_hdr_t *h; bam1_t *b;
    int seq_len; int64_t low = atoi(argv[2]), upper = atoi(argv[3]);
    if ((in = sam_open(argv[1], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[1]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[1]);
    b = bam_init1(); 

    if ((out = sam_open_format("-", "wb", NULL)) == NULL) err_fatal_simple("Cannot open \"-\"\n");
    if (sam_hdr_write(out, h) != 0) err_fatal_simple("Error in writing SAM header\n"); //sam header

    while (sam_read1(in, h, b) >= 0) {
        seq_len = _bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
        if ((low != -1 && (int64_t)seq_len < low) || (upper != -1 && (int64_t)seq_len > upper))
            continue;
        if (sam_write1(out, h, b) < 0) err_fatal_simple("Error in writing SAM record\n");
    }
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in); sam_close(out);
    return 0;
}

int fxt_filter_name(int argc, char* argv[]) {
    int i, c, exact_match=0, sub_match=0, input_list=0; char name[1024], sub_name[1024];
    int name_n=0, name_m = 10;
    char **name_array = (char**)_err_malloc(name_m * sizeof(char*));
    for (i = 0; i < name_m; ++i)
        name_array[i] = (char*)_err_malloc(1024 * sizeof(char));
    // hash table for names
    khash_t(str) *h = kh_init(str); int absent;
        
    while ((c = getopt(argc, argv, "n:m:l")) >= 0) {
        switch (c) {
            case 'n': exact_match=1, strcpy(name, optarg); break;
            case 'm': sub_match=1, strcpy(sub_name, optarg); break;
            case 'l': input_list = 1; break;
            default: err_printf("Error, unknown option: -%c %s\n", c, optarg);
        }
    }
    if (exact_match + sub_match != 1 || argc - optind != 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools filter-name [-n name] [-m sub-name] [-l] <in.fa/fq> > <out.fa/fq>\n");
        fprintf(stderr, "      -n [STR]    only output read with specified name.\n");
        fprintf(stderr, "      -m [STR]    only output read whose name or comment contain specified string.\n");
        fprintf(stderr, "      -l          input a list of names or sub-names with a list file, each line is a name or sub-name. [False]\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    char *fn = argv[optind];
    int is_fq = is_fastq(fn);
    if (is_fq < 0) err_fatal_simple("Unexpected input file format.");

    if (input_list) {
        FILE *fp; char line[1024];
        // read name/sub_name
        if (exact_match) fp = xopen(name, "r");
        else fp = xopen(sub_name, "r");
        while (fgets(line, 1024, fp) != NULL) {
            if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
            if (name_n == name_m) {
                name_m = name_n+1; kroundup32(name_m);
                name_array = (char**)_err_realloc(name_array, name_m * sizeof(char*));
                for (i = name_n; i < name_m; ++i)
                    name_array[i] = (char*)_err_malloc(1024 * sizeof(char));
            }
            strcpy(name_array[name_n], strtok(line, " \t"));

            khint_t pos = kh_put(str, h, name_array[name_n], &absent);
            if (absent) kh_val(h, pos) = name_n++;
        }
        err_fclose(fp);
    } else { // only 1 name
        if (exact_match) strcpy(name_array[name_n], name);
        else strcpy(name_array[name_n], sub_name);
        khint_t pos = kh_put(str, h, name_array[name_n], &absent);
        if (absent) kh_val(h, pos) = name_n++;
    }
    gzFile infp = xzopen(fn, "r"); FILE *out = stdout; kseq_t *seq = kseq_init(infp);
    while (kseq_read(seq) >= 0) {
        int hit = 0;
        if (exact_match) { // exact match
            khint_t pos = kh_get(str, h, seq->name.s);
            if (pos == kh_end(h)) hit = 0;
            else hit = 1;
        } else { // sub match
            if (seq->comment.l > 0) {
                for (i = 0; i < name_n; ++i) {
                    if (strstr(seq->name.s, name_array[i]) != NULL || strstr(seq->comment.s, name_array[i]) != NULL) {
                        hit = 1; break;
                    }
                }
            } else {
                for (i = 0; i < name_n; ++i) {
                    if (strstr(seq->name.s, name_array[i]) != NULL) {
                        hit = 1; break;
                    }
                }
            }
        }
        if (hit) print_seq(out, seq); 
    }
    err_fclose(out); kseq_destroy(seq); err_gzclose(infp);
    kh_destroy(str, h);
    for (i = 0; i < name_m; ++i) free(name_array[i]); free(name_array);
    return 0;
}

int fxt_filter_qual(int argc, char* argv[]) {
    int q_offset=33, c; 
    while ((c = getopt(argc, argv, "6")) >= 0) {
        switch (c) {
            case '6': q_offset=64; break;
            default: err_printf("Error, unknown option: -%c %s\n", c, optarg);
        }
    }
    if (argc - optind != 2) 
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools filter-qual <in.fa/fq> <min_qual> > <out.fa/fq>\n");
        fprintf(stderr, "      -6          quality value offset is 64. Default is 33. [False]\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[optind], "r"); float min_q = atof(argv[optind+1]);
    FILE *out = stdout;
    kseq_t *seq = kseq_init(infp);

    size_t i;
    while (kseq_read(seq) >= 0) {
        float sum_q = 0;
        for (i = 0; i < seq->qual.l; ++i) {
            sum_q += (int)seq->qual.s[i] - q_offset;
        }
        if (sum_q >= min_q * seq->qual.l) print_seq(out, seq); 
    }

    err_fclose(out);
    kseq_destroy(seq);
    err_gzclose(infp);
    return 0;
}

int fxt_filter_bam_name(int argc, char *argv[]) {
    int c, exact_match=0, sub_match=0, input_list=0; char name[1024], sub_name[1024];
    while ((c = getopt(argc, argv, "n:m:l")) >= 0) {
        switch (c) {
            case 'n': exact_match=1, strcpy(name, optarg); break;
            case 'm': sub_match=1, strcpy(sub_name, optarg); break;
            case 'l': input_list=1; break;
            default: err_printf("Error, unknown option: -%c %s\n", c, optarg);
        }
    }
    if (exact_match + sub_match != 1 || argc - optind != 1) 
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools filter-bam-name [-n name] [-m sub-name] [-l] <in.bam/sam> > <out.bam>\n");
        fprintf(stderr, "      -n [STR]    only output bam record with specified read name.\n");
        fprintf(stderr, "      -m [STR]    only output bam record whose read name contain specified string.\n");
        fprintf(stderr, "      -l          input a list of names or sub-names with a list file, each line is a name or sub-name. [False]\n");
        fprintf(stderr, "                  For example, \'fxtools fbn -n name.list in.bam > out.bam\'\n");
        fprintf(stderr, "\n");
        exit(-1);
    }

    int i, name_n=0, name_m = 10;
    char **name_array = (char**)_err_malloc(name_m * sizeof(char*));
    for (i = 0; i < name_m; ++i) 
        name_array[i] = (char*)_err_malloc(1024 * sizeof(char));
    // for fast exact match
    khash_t(str) *hash = kh_init(str); int absent;

    samFile *in, *out; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1(); 

    if ((out = sam_open_format("-", "wb", NULL)) == NULL) err_fatal_simple("Cannot open \"-\"\n");
    if (sam_hdr_write(out, h) != 0) err_fatal_simple("Error in writing SAM header\n"); //sam header

    char qname[1024];
    if (input_list) {
        // read name/sub_name
        FILE *fp;
        if (exact_match) {
            fp = xopen(name, "r");
        } else {
            fp = xopen(sub_name, "r");
        }
        char line[1024];
        while (fgets(line, 1024, fp) != NULL) {
            if (line[strlen(line)-1] == '\n')
                line[strlen(line)-1] = '\0';
            if (name_n == name_m) {
                name_m = name_n+1; kroundup32(name_m);
                name_array = (char**)_err_realloc(name_array, name_m * sizeof(char*));
                for (i = name_n; i < name_m ; ++i)
                    name_array[i] = (char*)_err_malloc(1024 * sizeof(char));
            }
            strcpy(name_array[name_n], line);
            khint_t pos = kh_put(str, hash, name_array[name_n], &absent);
            if (absent) kh_val(hash, pos) = name_n++;
        }
        err_fclose(fp);
    } else {
        if (exact_match) strcpy(name_array[name_n], name);
        else strcpy(name_array[name_n], sub_name);
        khint_t pos = kh_put(str, hash, name_array[name_n], &absent);
        if (absent) kh_val(hash, pos) = name_n++;
    }
    while (sam_read1(in, h, b) >= 0) {
        int hit = 0;
        strcpy(qname, bam_get_qname(b));
        if (exact_match) {
            khint_t pos = kh_get(str, hash, qname);
            if (pos == kh_end(hash)) hit = 0;
            else hit = 1;
        } else {
            for (i = 0; i < name_n; ++i) {
                if (strstr(qname, name_array[i]) != NULL) {
                    hit = 1;
                    break;
                }
            }
        }
        if (hit) {
            if (sam_write1(out, h, b) < 0) err_fatal_simple("Error in writing SAM record\n");
        }
    }
    for(i = 0; i < name_m; ++i) free(name_array[i]); free(name_array);
    kh_destroy(str, hash);
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in); sam_close(out);
    return 0;
}

int fxt_split_fx(int argc, char *argv[]) {
    int i;
    if (argc != 4) {
        err_printf("\nUsage: fxtools split-fx <in.fa/q> <N> <out_dir>\n\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[1], "r"); int n_files = atoi(argv[2]); char *out_dir = strdup(argv[3]);
    kseq_t *seq = kseq_init(infp);
    // file names of split files
    char *base = basename(argv[1]), out_name[1024];

    FILE **outfp = (FILE**)_err_malloc(n_files * sizeof(FILE*));
    for (i = 0; i < n_files; ++i) {
        sprintf(out_name, "%s/%s.%d", out_dir, base, i+1);
        outfp[i] = fopen(out_name, "w");
    }
    int read_i = 0; FILE *fp;
    while (kseq_read(seq) >= 0) {
        fp = outfp[read_i % n_files];
        print_seq(fp, seq);
        // fprintf(fp, ">%s\n", seq->name.s);
        // fprintf(fp, "%s\n", seq->seq.s);
        read_i++;
    }

    for (i = 0; i < n_files; ++i) {
        fclose(outfp[i]);
    } free(outfp);
    err_gzclose(infp);
    kseq_destroy(seq);
    return 0;
}

int fxt_fq2fa(int argc, char *argv[]) {
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools fq2fa <in.fq> > <out.fa>\n\n");
        exit(-1);
    } 
    gzFile infp = xzopen(argv[1], "r");
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = stdout;

    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, ">%s", seq->name.s);
        if (seq->comment.l > 0) fprintf(outfp, " %s", seq->comment.s);
        fprintf(outfp, "\n");

        fprintf(outfp, "%s\n", seq->seq.s);
    }

    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

int fxt_fa2fq(int argc, char *argv[]) {
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools fa2fq <in.fa> > <out.fq>\n\n");
        exit(-1);
    } 
    gzFile infp = xzopen(argv[1], "r");

    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = stdout;

    int64_t i;
    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, "@%s", seq->name.s);
        if (seq->comment.l > 0) fprintf(outfp, " %s", seq->comment.s);
        fprintf(outfp, "\n");

        fprintf(outfp, "%s\n", seq->seq.s);
        fprintf(outfp, "+\n");
        for (i = 0; i < (int64_t)seq->seq.l; ++i) fprintf(outfp, "I");
        fprintf(outfp, "\n");
    }

    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

int fxt_re_co(int argc, char *argv[]) {
    if (argc != 2)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools re-co in.fa/fq > out.fa\n"); fprintf(stderr, "\n");
        return 1;
    }
    gzFile readfp;
    kseq_t *read_seq;
    int seq_len = 1024;
    char *seq = (char*)malloc(seq_len*sizeof(char));
    int8_t *seq_n = (int8_t*)malloc(seq_len*sizeof(int8_t));
    FILE *out = stdout;

    readfp = xzopen(argv[1], "r");
    read_seq = kseq_init(readfp);

    int len, i;

    while ((len = kseq_read(read_seq)) > 0)
    {
       if (len > seq_len)
       {
           seq_len = len; kroundup32(seq_len);
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
       if (read_seq->qual.l > 0) {
           fprintf(out, "@%s_reverse_complementary", read_seq->name.s);
           if (read_seq->comment.l > 0) fprintf(out, " %s", read_seq->comment.s);
           fprintf(out, "\n");
           for (i = len - 1; i>=0; i--)
           {
               if (seq_n[i] != 4) 
                   fprintf(out, "%c", nt_char[3-(int)seq_n[i]]);
                else
                    fprintf(out, "N");
           }
           fprintf(out, "\n+\n");
           for (i = len - 1; i>=0; i--) fprintf(out, "%c", read_seq->qual.s[i]);
           fprintf(out, "\n");

       } else {
           fprintf(out, ">%s_reverse_complementary", read_seq->name.s);
           if (read_seq->comment.l > 0) fprintf(out, " %s", read_seq->comment.s);
           fprintf(out, "\n");
           for (i = len - 1; i>=0; i--)
           {
               if (seq_n[i] != 4) 
                   fprintf(out, "%c", nt_char[3-(int)seq_n[i]]);
                else
                    fprintf(out, "N");
           }
           fprintf(out, "\n");
       }
    }

    err_gzclose(readfp);
    kseq_destroy(read_seq);
    err_fclose(out);

    return 0;
}

int fxt_dna2rna(int argc, char *argv[]) {
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools dna2rna <in.fa/fq> > <out.fa/fq>\n\n");
        exit(-1);
    } 
    gzFile infp = xzopen(argv[1], "r");
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = stdout;

    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, ">%s", seq->name.s);
        if (seq->comment.l > 0) fprintf(outfp, " %s", seq->comment.s);
        fprintf(outfp, "\n");

        size_t i;
        for (i = 0; i < seq->seq.l; ++i) {
            if (seq->seq.s[i] == 'T') fprintf(outfp, "U");
            else fprintf(outfp, "%c", seq->seq.s[i]);
        }
        fprintf(outfp, "\n");
    }

    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

int fxt_rna2dna(int argc, char *argv[]) {
    if (argc != 2)
    {
        fprintf(stderr, "\n"); fprintf(stderr, "Usage: fxtools rna2dna <in.fa/fq> > <out.fa/fq>\n\n");
        exit(-1);
    } 
    gzFile infp = xzopen(argv[1], "r");
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = stdout;

    while (kseq_read(seq) >= 0)
    {
        fprintf(outfp, ">%s", seq->name.s);
        if (seq->comment.l > 0) fprintf(outfp, " %s", seq->comment.s);
        fprintf(outfp, "\n");

        size_t i;
        for (i = 0; i < seq->seq.l; ++i) {
            if (seq->seq.s[i] == 'U') fprintf(outfp, "T");
            else fprintf(outfp, "%c", seq->seq.s[i]);
        }
        fprintf(outfp, "\n");
    }

    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

int fxt_seq_dis(int argc, char *argv[]) {
    if (argc != 5)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools seq-display <in.fa/fq> <chr/read_name> <start_pos(1-based)> <end_pos>\n");
        fprintf(stderr, "       use negative coordinate to indicate later part of sequence. (e.g., -1 for last bp)\n");
        fprintf(stderr, "\n"); 
        exit(-1);
    }
    int is_fq = is_fastq(argv[1]);
    if (is_fq < 0) err_fatal_simple("Unexpected input file format.");

    char *fai_fn = (char*)_err_malloc(strlen(argv[1])+5);
    strcpy(fai_fn, argv[1]); strcat(fai_fn, ".fai");
    faidx_t *fai = fai_load3_format(argv[1], fai_fn, NULL, FAI_CREATE, is_fq == 0 ? FAI_FASTA : FAI_FASTQ);
    fprintf(stderr, "Loading/building fai index %s.fai ... ", argv[1]);
    if ( !fai ) {
        fprintf(stderr, "\nCould not load/build fai index %s.fai\n", argv[1]);
        return EXIT_FAILURE;
    }
    fprintf(stderr, "done!\n");

    int exit_status = EXIT_SUCCESS;

    char reg[1024], chr[102]; int start, end;
    strcpy(chr, argv[2]);
    start = atoi(argv[3]); end = atoi(argv[4]);
    int tot_len = faidx_seq_len(fai, chr);
    if (start < 0) start = tot_len + start + 1;
    if (end < 0) end = tot_len + end + 1;
    sprintf(reg, "%s:%d-%d", chr, start, end);
    printf(">%s\n", reg);
    int seq_len;
    char *seq = fai_fetch(fai, reg, &seq_len);
    if ( seq_len < 0 ) {
        err_printf("Failed to fetch sequence in %s\n", reg);
        exit_status = EXIT_FAILURE;
        return exit_status;
    }
    size_t i, seq_sz = seq_len;
    for (i=0; i<seq_sz; i+=60)
    {
        size_t len = i + 60 < seq_sz ? 60 : seq_sz - i;
        if (fwrite(seq + i, 1, len, stdout) < len ||
                putchar('\n') == EOF) {
            err_fatal_simple("failed to write output");
        }
    }
    free(seq);
    fai_destroy(fai);
    return exit_status;
}

int fxt_cigar_parse(int argc, char *argv[]) {
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
    setlocale(LC_NUMERIC, "");
    fprintf(stdout, "\nCigar length:\n");
    for (i = 0; i < 11; ++i)
    {
        if (op[i] != 0) fprintf(stdout, "%'16ld %c\n", op[i], CIGAR_STR[i]);
    }
    fprintf(stdout, "\nseq-len: %'7d\nref-len: %'7d\n", seq_len, ref_len);
    return 0;
}

int int_cmp(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

void print_len_stats(char *fn, int *len, int n) {
    setlocale(LC_NUMERIC, "");
    int i, min_len=INT32_MAX, max_len=INT32_MIN, n50_len=0, mean_len = 0;
    long long tot_len = 0, n50_tot_len = 0;

    qsort(len, n, sizeof(int), int_cmp);
    for (i = 0; i < n; ++i) {
        tot_len += len[i];
        if (len[i] > max_len) max_len = len[i];
        if (len[i] < min_len) min_len = len[i];
    }
    if (n == 0) mean_len = 0;
    else mean_len = tot_len / n;

    for (i = n-1; i >= 0; --i) {
        n50_tot_len += len[i];
        if (n50_tot_len >= tot_len / 2) {
            n50_len = len[i];
            break;
        }
    }
    fprintf(stderr, "== \'%s\' read length stats ==\n", fn);
    fprintf(stderr, "Total reads\t%'16d\n", n);
    fprintf(stderr, "Total bases\t%'16lld\n", tot_len);
    fprintf(stderr, "Mean length\t%'16d\n", mean_len);
    fprintf(stderr, "Min. length\t%'16d\n", min_len);
    fprintf(stderr, "Max. length\t%'16d\n", max_len);
    fprintf(stderr, "N-50 length\t%'16d\n", n50_len);
}

int fxt_len_parse(int argc, char *argv[]) {
    if (argc < 2)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools length-parse <in.fa/fq/len>\n");
        fprintf(stderr, "\n"); 
        exit(-1);
    }
    int i, *len, n, m;
    n = 0, m = 1000;
    len = (int*)_err_malloc(m * sizeof(int));
    for(i = 1; i < argc; ++i) {
        int not_len_file = 1;
        n = 0;
        gzFile infp = xzopen(argv[i], "r");
        char buff[1024]; int buff_n;
        // check if input file is .len file
        if (strcmp(argv[i], "-") != 0) {
            buff_n = gzread(infp, buff, 10);
            if (buff[0] != '>' && buff[0] != '@')
                not_len_file = 0;
            err_gzclose(infp);
            infp  = xzopen(argv[i], "r");
        } else
            not_len_file = 1;
        if (not_len_file) {
            kseq_t *seq;
            seq = kseq_init(infp);
            while (kseq_read(seq) >= 0)
            {
                fprintf(stdout, "%s\t%d\n", seq->name.s, (int)seq->seq.l);       
                if (n == m) {
                    m = n+1; kroundup32(m);
                    len = (int*)_err_realloc(len, m * sizeof(int));
                }
                len[n++] = seq->seq.l;
            }
            kseq_destroy(seq);
        } else { // is .len file
            int seq_len;
            while (gzgets(infp, buff, 1024) != NULL) {
                sscanf(buff, "%*s %d", &seq_len);
                if (n == m) {
                    m = n+1; kroundup32(m);
                    len = (int*)_err_realloc(len, m * sizeof(int));
                }
                len[n++] = seq_len; 
            }
        }
        err_gzclose(infp);
        print_len_stats(argv[i], len, n);
    }
    free(len);
    return 0;
}

int fxt_merge_filter_fa(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools merge-fil-fa <in.fa> [N] > <out.fa/fq>\n");
        fprintf(stderr, "         optional: use N to separate merged sequences\n");
        fprintf(stderr, "         only work with fasta file.\n"); 
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[1], "r");
    kseq_t *seq = kseq_init(infp);
    FILE *outfp = stdout;
    char read_name[1024]; 
    char *read_seq = (char*)calloc(10, 1);
    int w_seq_n=0;

    int i, j, i1=-1, i2=-1, n=0;
    int *len=(int*)_err_malloc(100 * sizeof(int)), *tmp=(int*)_err_malloc(100 * sizeof(int));
    char **name=(char**)_err_malloc(100 * sizeof(char*));
    for (i = 0; i < 100; i++) name[i] = (char*)_err_malloc(20 * sizeof(char));
    char sep[5]; 
    if (argc == 3) strcpy(sep, "N");
    else strcpy(sep, "");
    while (kseq_read(seq) >= 0)
    {
        if (strcmp(seq->name.s, read_name) == 0) {
            w_seq_n += seq->seq.l;
            strcpy(name[n], seq->comment.s);
            len[n++] = seq->seq.l;
            read_seq = (char*)realloc(read_seq, w_seq_n+seq->seq.l);
            strcat(read_seq, seq->seq.s);
            read_seq[w_seq_n] = 0;
        } else {
            if (w_seq_n > 0) {
                // cal i1 and i2
                i1 = -1, i2 = -1;
                for (i = 0; i < n; ++i) tmp[i] = len[i];
                qsort(tmp, n, sizeof(int), int_cmp);
                if (n > 1 && n <= 10) {
                    fprintf(outfp, ">%s r1", read_name);
                    for (i = 0; i < n; ++i) {
                        if (len[i] == tmp[(n-1)/2]) {
                            i1 = i;
                            fprintf(outfp, " %s", name[i]);
                            break;
                        }
                    }
                    fprintf(outfp, "\n");
                } else if (n > 10) {
                    fprintf(outfp, ">%s r2", read_name);
                    for (i = 0; i < n; ++i) {
                        if (len[i] == tmp[(n-1)/3]) {
                            i1 = i;
                            fprintf(outfp, " %s", name[i]);
                        } else if (len[i] == tmp[(n-1)*2/3]){
                            i2 = i;
                            fprintf(outfp, " %s", name[i]);
                        }
                    }
                    fprintf(outfp, "\n");
                } else {
                    fprintf(outfp, ">%s\n", read_name);
                }
                // filter with i1 and i2
                int start = 0, end = 0, first = 0;
                for (i = 0; i < n; ++i) {
                    end += len[i];
                    if (i != i1 && i != i2) {
                        if (first) fprintf(outfp, "%s", sep);
                        for (j = start; j < end; ++j)
                            fprintf(outfp, "%c", read_seq[j]);
                        first = 1;
                    }
                    start = end;

                }
                fprintf(outfp, "\n");
            }
            n = 0;
            w_seq_n = seq->seq.l;
            strcpy(name[n], seq->comment.s);
            len[n++] = seq->seq.l;
            strcpy(read_name, seq->name.s);
            read_seq = (char*)realloc(read_seq, w_seq_n+seq->seq.l);
            strcpy(read_seq, seq->seq.s);
            read_seq[w_seq_n] = 0;
        }
    }
    if (w_seq_n > 0) { // last read
        // cal i1 and i2
        i1 = -1, i2 = -1;
        for (i = 0; i < n; ++i) tmp[i] = len[i];
        qsort(tmp, n, sizeof(int), int_cmp);
        if (n > 1 && n < 10) {
            fprintf(outfp, ">%s r1", read_name);
            for (i = 0; i < n; ++i) {
                if (len[i] == tmp[(n-1)/2]) {
                    i1 = i;
                    fprintf(outfp, " %s", name[i]);
                    break;
                }
            }
            fprintf(outfp, "\n");
        } else if (n > 10){
            fprintf(outfp, ">%s r2", read_name);
            for (i = 0; i < n; ++i) {
                if (len[i] == tmp[(n-1)/3]) {
                    i1 = i;
                    fprintf(outfp, " %s", name[i]);
                } else if (len[i] == tmp[(n-1)*2/3]){
                    i2 = i;
                    fprintf(outfp, " %s", name[i]);
                }
            }
            fprintf(outfp, "\n");
        } else {
            fprintf(outfp, ">%s\n", read_name);
        }
        int start = 0, end = 0, first = 0;
        for (i = 0; i < n; ++i) {
            end += len[i];
            if (i != i1 && i != i2) {
                if (first) fprintf(outfp, "%s", sep);
                for (j = start; j < end; ++j)
                    fprintf(outfp, "%c", read_seq[j]);
                first = 1;
            }
            start = end;
        }
        fprintf(outfp, "\n");
    }
    kseq_destroy(seq);
    free(name); free(len); free(tmp); err_gzclose(infp); err_fclose(outfp);
    return 0;
}

int fxt_merge_fa(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools merge-fa <in.fa/fq> [N] > <out.fa/fq>\n");
        fprintf(stderr, "         optional: use N to separate merged sequences\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[1], "r");
    kseq_t *seq = kseq_init(infp);
    FILE *outfp = stdout;
    char read_name[1024]; 
    char *read_seq = (char*)calloc(10, 1);
    char *read_qual = (char*)calloc(10, 1);
    int w_seq_n=0;

    char sep[5];
    if (argc == 3) strcpy(sep, "N");
    else strcpy(sep, "");
    while (kseq_read(seq) >= 0)
    {
        if (strcmp(seq->name.s, read_name) == 0) {
            w_seq_n += seq->seq.l;
            read_seq = (char*)realloc(read_seq, w_seq_n+seq->seq.l);
            strcat(read_seq, sep);
            strcat(read_seq, seq->seq.s);
            read_seq[w_seq_n] = 0;
            if (seq->qual.l > 0) {
                read_qual = (char*)realloc(read_qual, w_seq_n+seq->seq.l);
                strcat(read_qual, "!");
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
    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

int fxt_duplicate_read(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools duplicate-read <in.fa/fq> <copy_number> > out.fa/fq\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[1], "r"); float copy_n = atof(argv[2]);
    kseq_t *seq = kseq_init(infp);
    FILE *outfp = stdout;
    size_t i; int j;

    while (kseq_read(seq) >= 0)
    {
        if (seq->qual.l > 0) { // fastq
            for (j = 0; j < copy_n; ++j) {
                fprintf(outfp, "@%s_%d", seq->name.s, j);
                if (seq->comment.l > 0) fprintf(outfp, " %s\n", seq->comment.s);
                else fprintf(outfp, "\n");
                for (i = 0; i < seq->seq.l; ++i) {
                    fprintf(outfp, "%c", seq->seq.s[i]);
                } fprintf(outfp, "\n");
                fprintf(outfp, "+\n");
                for (i = 0; i < seq->qual.l; ++i) {
                    fprintf(outfp, "%c", seq->qual.s[i]);
                } fprintf(outfp, "\n");
            }
        } else { // fasta
            for (j = 0; j < copy_n; ++j) {
                fprintf(outfp, ">%s_%d", seq->name.s, j);
                if (seq->comment.l > 0) fprintf(outfp, " %s\n", seq->comment.s);
                else fprintf(outfp, "\n");
                for (i = 0; i < seq->seq.l; ++i) {
                    fprintf(outfp, "%c", seq->seq.s[i]);
                } fprintf(outfp, "\n");
            }
        }
    }
       
    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

int fxt_duplicate_seq(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fxtools duplicate-seq <in.fa/fq> <copy_number> > out.fa/fq\n");
        fprintf(stderr, "\n");
        exit(-1);
    }
    gzFile infp = xzopen(argv[1], "r"); float copy_n = atof(argv[2]);
    kseq_t *seq = kseq_init(infp);
    FILE *outfp = stdout;
    size_t i;

    while (kseq_read(seq) >= 0)
    {
        if (seq->qual.l > 0) { // fastq
            fprintf(outfp, "@%s_copy:%.2f", seq->name.s, copy_n);
            if (seq->comment.l > 0) fprintf(outfp, " %s\n", seq->comment.s);
            else fprintf(outfp, "\n");
            for (i = 0; i < seq->seq.l * copy_n; ++i) {
                fprintf(outfp, "%c", seq->seq.s[i % seq->seq.l]);
            } fprintf(outfp, "\n");
            fprintf(outfp, "+\n");
            for (i = 0; i < seq->qual.l * copy_n; ++i) {
                fprintf(outfp, "%c", seq->qual.s[i % seq->qual.l]);
            } fprintf(outfp, "\n");
        } else { // fasta
            fprintf(outfp, ">%s_copy:%.2f", seq->name.s, copy_n);
            if (seq->comment.l > 0) fprintf(outfp, " %s\n", seq->comment.s);
            else fprintf(outfp, "\n");
            for (i = 0; i < seq->seq.l * copy_n; ++i) {
                fprintf(outfp, "%c", seq->seq.s[i % seq->seq.l]);
            } fprintf(outfp, "\n");
        }
    }
       
    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

int fxt_error_parse(int argc, char *argv[]) {
    setlocale(LC_NUMERIC, "");
    int c, parse_non_primary = 0;
    while ((c = getopt(argc, argv, "s")) >= 0) {
        switch (c) {
            case 's': parse_non_primary=1; break;
            default: err_printf("Error, unknown option: -%c %s\n", c, optarg);
        }
    }
    if (argc - optind != 1)
    {
        fprintf(stderr, "\n"); 
        fprintf(stderr, "Usage: fxtools error-parse <input.SAM/BAM/PAF/GAF> [-s] > error.out\n");
        fprintf(stderr, "         -s    include non-primary records in the output.\n\n");
        return 1;
    }
    fprintf(stdout, "#READ_NAME\tREAD_LEN\tUNMAP\tINS\tDEL\tMIS\tMATCH\tCLIP\tSKIP\tERR_RATE\n");
    char *qname;
    long long tol_n=0, unmap=0, is_primary=0, tol_len=0, tol_ins=0, tol_del=0, tol_mis=0, tol_match=0, tol_clip=0, tol_skip=0;
    int i, seq_len, unmap_flag=0, md, ins, del, mis, match, clip, skip;
    int equal, diff;

    if (check_suf(argv[optind], ".sam") || check_suf(argv[optind], ".bam")) {
        samFile *in; bam_hdr_t *h; bam1_t *b;
        if ((in = sam_open(argv[optind], "r")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
        b = bam_init1(); 

        while (sam_read1(in, h, b) >= 0) {
            if (!parse_non_primary && (b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY)) continue;
            tol_n++;
            unmap_flag = 0;
            md = 0, ins = 0, del = 0, mis = 0, match = 0, clip = 0, skip = 0;
            equal = 0, diff = 0;
            if (!bam_unmap(b)) {
                uint32_t *cigar = bam_get_cigar(b); int cigar_len = b->core.n_cigar;
                cigar_parse(b, cigar, cigar_len, &match, &mis, &ins, &del, &skip, &clip);
                seq_len = _bam_cigar2qlen(cigar_len, cigar);
            } else {
                unmap++;
                unmap_flag = 1;
                seq_len = b->core.l_qseq;
            }
            tol_len += seq_len; tol_ins += ins; tol_del += del; tol_mis += mis; tol_match += match; tol_clip += clip; tol_skip += skip;

            fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f%%\n", bam_get_qname(b), seq_len, unmap_flag, ins, del, mis, match, clip, skip, (ins+mis+del+0.1)/(match+ins+mis)*100);
        }
        bam_destroy1(b); sam_close(in); bam_hdr_destroy(h);
    } else if (check_suf(argv[optind], ".gaf") || check_suf(argv[optind], ".gaf.gz")) {
        gzFile gaf_fp = err_xzopen_core(__func__, argv[optind], "r"); kstream_t *ks = ks_init(gaf_fp);
        kstring_t line = {0,0,0}; int dret;
        while (ks_getuntil(ks, KS_SEP_LINE, &line, &dret) >= 0) {
            unmap_flag = 0; is_primary = 0;
            md = 0, ins = 0, del = 0, mis = 0, match = 0, clip = 0, skip = 0;
            // GAF: 
            //  0:qname, 1:qlen, 2:qs, 3: qe, 4:+/-, 5:path, 6:plen, 7:ps, 8:pe, 9:n_match, 10:n_aln, 11:mapq
            int l_aux, m_aux=0; uint8_t *aux=0, *info;
            char *deli_s, *info_s, *aux_s; int is_ok = 0;
            for (i = 0, deli_s = info_s = line.s;; ++deli_s) {
                if (*deli_s == 0 || *deli_s == '\t') {
                    int c = *deli_s;
                    *deli_s = 0;
                    if (i == 0) {
                        qname = info_s;
                    } else if (i == 1) {
                        seq_len = _strtol10(info_s, NULL);
                    } else if (i == 5) { // path
                        if (info_s[0] == '*') { // unmaped
                            unmap++;
                            unmap_flag = 1;
                            break;
                        }
                    } else if (i == 11) {
                        is_ok = 1; aux_s = deli_s + 1;
                        break;
                    }
                    if (c == 0) break;
                    ++i, info_s = deli_s + 1;
                }
            }
            if (is_ok) {
                l_aux = aux_parse(aux_s, &aux, &m_aux);
                info = aux_get(l_aux, aux, "tp");
                if (info == 0 || info[0] != 'A') err_fatal(__func__, "Error: no \"tp\" tag in line\n%s", line.s);
                if (!parse_non_primary && info[1] != 'P') continue;
                l_aux = aux_del(l_aux, aux, info);

                info = aux_get(l_aux, aux, "cg");
                if (info == 0 || info[0] != 'Z') err_fatal(__func__, "Error: no \"cg\" tag in line\n%s", line.s);
                char *cigar_str = (char*)info+1;
                cigar_str_parse(cigar_str, &match, &mis, &ins, &del, &skip, &clip);
                l_aux = aux_del(l_aux, aux, info);
                if (aux) free(aux);
            }
            fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f%%\n", qname, seq_len, unmap_flag, ins, del, mis, match, clip, skip, (0.0+ins+del+mis)/(match+ins+mis)*100);
            tol_len += seq_len; tol_ins += ins; tol_del += del; tol_mis += mis; tol_match += match; tol_clip += clip; tol_skip += skip;
            tol_n++;
        }
        ks_destroy(ks); err_gzclose(gaf_fp); if (line.m) free(line.s);
    } else if (check_suf(argv[optind], ".paf") || check_suf(argv[optind], ".paf.gz")) {
        gzFile paf_fp = err_xzopen_core(__func__, argv[optind], "r"); kstream_t *ks = ks_init(paf_fp);
        kstring_t line = {0,0,0}; int dret;
        while (ks_getuntil(ks, KS_SEP_LINE, &line, &dret) >= 0) {
            unmap_flag = 0; is_primary = 0;
            md = 0, ins = 0, del = 0, mis = 0, match = 0, clip = 0, skip = 0;
            // PAF: 
            //  0:qname, 1:qlen, 2:qs(0), 3: qe(0), 4:+/-, 5:tname, 6:tlen, 7:ts, 8:te, 9:n_match, 10:n_aln, 11:mapq
            int l_aux, m_aux=0; uint8_t *aux=0, *info;
            char *deli_s, *info_s, *aux_s; int is_ok = 0;
            for (i = 0, deli_s = info_s = line.s;; ++deli_s) {
                if (*deli_s == 0 || *deli_s == '\t') {
                    int c = *deli_s;
                    *deli_s = 0;
                    if (i == 0) {
                        qname = info_s;
                    } else if (i == 1) {
                        seq_len = _strtol10(info_s, NULL);
                    } else if (i == 5) { // path
                        if (info_s[0] == '*') { // unmaped
                            unmap++;
                            unmap_flag = 1;
                            break;
                        }
                    } else if (i == 11) {
                        is_ok = 1; aux_s = deli_s + 1;
                        break;
                    }
                    if (c == 0) break;
                    ++i, info_s = deli_s + 1;
                }
            }
            if (is_ok) {
                l_aux = aux_parse(aux_s, &aux, &m_aux);
                info = aux_get(l_aux, aux, "tp");
                if (info == 0 || info[0] != 'A') err_fatal(__func__, "Error: no \"tp\" tag in line\n%s", line.s);
                if (!parse_non_primary && info[1] != 'P') continue;
                l_aux = aux_del(l_aux, aux, info);

                int32_t NM = -1;
                if ((info = aux_get(l_aux, aux, "NM"))) {
                    if (info[0] != 'i') err_fatal(__func__, "Error in \"NM\" tag in line\n%s", line.s);
                    NM = *(int32_t*)(info+1);
                }

                info = aux_get(l_aux, aux, "cg");
                if (info == 0 || info[0] != 'Z') err_fatal(__func__, "Error: no \"cg\" tag in line\n%s", line.s);
                char *cigar_str = (char*)info+1;
                cigar_str_parse(cigar_str, &match, &mis, &ins, &del, &skip, &clip);
                if (ins+del+mis != NM) {
                    mis = NM - ins -del;
                    match -= mis;
                }
                l_aux = aux_del(l_aux, aux, info);
                if (aux) free(aux);
            }
            fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f%%\n", qname, seq_len, unmap_flag, ins, del, mis, match, clip, skip, (0.0+ins+del+mis)/(match+ins+mis)*100);
            tol_len += seq_len; tol_ins += ins; tol_del += del; tol_mis += mis; tol_match += match; tol_clip += clip; tol_skip += skip;
            tol_n++;
        }
        ks_destroy(ks); err_gzclose(paf_fp); if (line.m) free(line.s);
    } else {
        err_fatal(__func__, "Unexpected file format: %s\n", argv[optind]);
    }

    fprintf(stdout, "%s\t%'lld\t%'lld\t%'lld\t%'lld\t%'lld\t%'lld\t%'lld\t%'lld\t%.1f%%\n", "#TOTAL", tol_len, unmap, tol_ins, tol_del, tol_mis, tol_match, tol_clip, tol_skip, (tol_ins+tol_del+tol_mis+0.0)/(tol_match+tol_ins+tol_mis)*100);
    fprintf(stdout, "#Total mapped read: %'lld (%.1f%%)\n#Total unmapped read: %'lld\n#Total read: %'lld\n#Error rate: %.1f%%\n", tol_n-unmap, (tol_n-unmap+0.0)/ tol_n * 100, unmap, tol_n, (tol_ins+tol_del+tol_mis+0.0)/(tol_match+tol_ins+tol_mis) * 100); // no tol_del
    return 0;
}

int get_read_pos(int pos, uint32_t *cigar, int n_cigar, int ref_pos/*1-based*/, int read_pos) {
    int i;
    for (i = 0; i < n_cigar; ++i) {
        int l = cigar[i]>>4, op=cigar[i]&0xf;
        if (op == BAM_CREF_SKIP || op == BAM_CDEL) {
            if (ref_pos + l > pos) {
                return read_pos;
            }
            ref_pos += l;
        } else if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // ref:1, read:1
            if (ref_pos + l > pos) {
                return read_pos + (pos-ref_pos);
            }
            read_pos += l, ref_pos += l;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // ref:0, read:1
            read_pos += l;
        }
    }
    return read_pos;
}

int bam2bed_core(bam1_t *b, bam_hdr_t *h) {
    int n_cigar = b->core.n_cigar; uint32_t *cigar = bam_get_cigar(b);
    int is_rev = bam_is_rev(b);
    int i, read_s=1, read_e = 0, ref_s = b->core.pos+1, ref_e = ref_s - 1, tmp_s, tmp_e;
    int qlen = _bam_cigar2qlen(n_cigar, cigar);
    for (i = 0; i < n_cigar; ++i) {
        int l = cigar[i]>>4, op=cigar[i]&0xf;
        if (op == BAM_CREF_SKIP) {
            if (is_rev) {
                tmp_e = qlen + 1 - read_s;
                tmp_s = qlen + 1 - read_e;
            } else {
                tmp_s = read_s;
                tmp_e = read_e;
            }
            fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\n", h->target_name[b->core.tid], ref_s, ref_e, bam_get_qname(b), 0, "+-"[is_rev], tmp_s, tmp_e);
            ref_s = ref_e + l + 1;
            ref_e += l;
            read_s = read_e + 1;
        } else if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // ref:1, read:1
            read_e += l, ref_e += l;
        } else if (op == BAM_CDEL) { // ref:1, read:0
            ref_e += l;
        } else if (op == BAM_CINS) { // ref:0, read:1
            read_e += l;
        } else if (op == BAM_CSOFT_CLIP) { // ref: 0, read: 1
            if (i == 0) {
                read_s += l;
                read_e += l;
            }
        }
    }
    if (is_rev) {
        tmp_e = qlen + 1 - read_s;
        tmp_s = qlen + 1 - read_e;
    } else {
        tmp_s = read_s;
        tmp_e = read_e;
    }
    fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\n", h->target_name[b->core.tid], ref_s, ref_e, bam_get_qname(b), 0, "+-"[is_rev], tmp_s, tmp_e);

    return 0;
}

int fxt_bam2bed(int argc, char *argv[]) {
    if (argc != 2) {
        err_printf("Usage: fxtools bam2bed in.bam > out.bed\n");
        err_printf("\n");
        return 0;
    }
    char bamfn[1024];
    strcpy(bamfn, argv[1]);
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(bamfn, "rb")) == NULL) err_fatal(__func__, "fail to open \"%s\"\n", bamfn);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "fail to read header for \"%s\"\n", bamfn);
    b = bam_init1();
    int r;
    while ((r = sam_read1(in, h, b)) >= 0) {
        if (b->core.flag & BAM_FUNMAP) continue;
        bam2bed_core(b, h);
    }
    bam_hdr_destroy(h); sam_close(in); bam_destroy1(b);
    return 0;
}

int fxt_bam2fx(int argc, char *argv[]) {
    int c, out_ref_base = 0, only_mapped=0, fq_out=0, include_non_primary=0;
    while ((c = getopt(argc, argv, "rsmq")) >= 0) {
        switch (c) {
            case 'r': out_ref_base = 1; break;
            case 's': include_non_primary = 1; break;
            case 'm': only_mapped=1; break;
            case 'q': fq_out= 1; break;
            default: err_printf("Error, unknown option: -%c %s\n", c, optarg); break;
        }
    }
    if (argc - optind < 1) {
        err_printf("Usage: fxtools bam2fx in.bam [-rsmq] > out.fa/fq\n");
        err_printf("         -r    output reference sequence base, instead of raw read sequence base.\n");
        err_printf("         -s    include non-primary records in the output.\n");
        err_printf("         -m    only output mapped sequence (discard clipped unmapped sequence).\n");
        err_printf("         -q    output fastq format.\n");
        err_printf("\n");
        return 0;
    }
    if (only_mapped==0 && include_non_primary==1) {
        fprintf(stderr, "Warning: For supplementary/secondary alignment records, only mapped seqeunce are output.\n");
    }
    char bamfn[1024];
    strcpy(bamfn, argv[optind]);
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(bamfn, "rb")) == NULL) err_fatal(__func__, "fail to open \"%s\"\n", bamfn);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "fail to read header for \"%s\"\n", bamfn);
    b = bam_init1();
    int i, r;
    uint8_t *bam_bseq; char *bam_qual;
    while ((r = sam_read1(in, h, b)) >= 0) {
        if (b->core.flag & BAM_FUNMAP) continue;
        if (include_non_primary == 0 && (b->core.flag & BAM_FSUPPLEMENTARY || b->core.flag & BAM_FSECONDARY)) continue;
        bam_bseq = bam_get_seq(b); bam_qual = bam_get_qual(b);
        int start, end, len = b->core.l_qseq, left_clip = 0, right_clip = 0;
        int is_rev = (b->core.flag & BAM_FREVERSE) != 0;
        if (only_mapped && !(b->core.flag & BAM_FSUPPLEMENTARY) && !(b->core.flag & BAM_FSECONDARY)) { // discard clipping sequence
            int n_cigar = b->core.n_cigar; uint32_t *cigar = bam_get_cigar(b);
            for (i = 0; i < n_cigar; ++i) {
                int op = bam_cigar_op(cigar[i]);
                if (i == 0 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
                    left_clip = bam_cigar_oplen(cigar[i]);
                }
                if (i == n_cigar-1 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
                    right_clip = bam_cigar_oplen(cigar[i]);
                }
            }
            start = left_clip; end = len - right_clip;
        } else {
            start = 0; end = len;
        }
        if (out_ref_base == 0 && is_rev == 1) { // rc seq
            if (fq_out) { // fastq
                fprintf(stdout, "@%s %d %d %d\n", bam_get_qname(b), len, start, end);
                for (i = end-1; i >= start; --i) {
                    fprintf(stdout, "%c", bam_rc_qseq_table[bam_seqi(bam_bseq, i)]);
                } fprintf(stdout, "\n");
                fprintf(stdout, "+\n");
                if (bam_qual[0] == 0xff) fprintf(stdout, '*');
                else for (i = end-1; i >= start; --i) fprintf(stdout, "%c", bam_qual[i] + 33);
                fprintf(stdout, "\n");
            } else { // fasta
                fprintf(stdout, ">%s %d %d %d\n", bam_get_qname(b), len, start, end);
                for (i = end-1; i >= start; --i) {
                    fprintf(stdout, "%c", bam_rc_qseq_table[bam_seqi(bam_bseq, i)]);
                } fprintf(stdout, "\n");
            }
        } else {
            if (fq_out) { // fastq
                fprintf(stdout, "@%s %d %d %d\n", bam_get_qname(b), len, start, end);
                for (i = start; i < end; ++i) {
                    fprintf(stdout, "%c", bam_qseq_table[bam_seqi(bam_bseq, i)]);
                } fprintf(stdout, "\n");
                fprintf(stdout, "+\n");
                if (bam_qual[0] == 0xff) fprintf(stdout, '*');
                else for (i = start; i < end; ++i) fprintf(stdout, "%c", bam_qual[i] + 33);
                fprintf(stdout, "\n");
            } else { // fasta
                fprintf(stdout, ">%s %d %d %d\n", bam_get_qname(b), len, start, end);
                for (i = start; i < end; ++i) {
                    fprintf(stdout, "%c", bam_qseq_table[bam_seqi(bam_bseq, i)]);
                } fprintf(stdout, "\n");
            }
        }

    }
    bam_hdr_destroy(h); sam_close(in); bam_destroy1(b); 
    return 0;
}

int fxt_sam_flag(int argc, char *argv[]) {
    if (argc != 3) {
        err_printf("Usage: fxtools sam-flag in.sam/bam flag > flag.out\n");
        err_printf("\n");
        return 0;
    }
    char bamfn[1024], flag[10];
    strcpy(bamfn, argv[1]); strcpy(flag, argv[2]);
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(bamfn, "rb")) == NULL) err_fatal(__func__, "fail to open \"%s\"\n", bamfn);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "fail to read header for \"%s\"\n", bamfn);
    b = bam_init1();
    int r;
    while ((r = sam_read1(in, h, b)) >= 0) {
        uint8_t *p = bam_aux_get(b, flag);
        if (p == 0) {
            fprintf(stderr, "%s has No \"%s\" tag.\n", bam_get_qname(b), flag);
        } else {
            if (*p == 'i' || *p == 'C' || *p == 'c' || *p == 's' || *p == 'S') {
                int i = bam_aux2i(p);
                fprintf(stdout, "%s\t%d\n", bam_get_qname(b), i);
            } else if (*p == 'd' || *p == 'f') {
                double d = bam_aux2f(p);
                fprintf(stdout, "%s\t%lf\n", bam_get_qname(b), d);
            } else if (*p == 'A') { 
                char c = bam_aux2A(p);
                fprintf(stdout, "%s\t%c\n", bam_get_qname(b), c);
            } else if (*p == 'Z' || *p == 'H') { 
                char *c = bam_aux2Z(p);
                fprintf(stdout, "%s\t%s\n", bam_get_qname(b), c);
            }
        }
    }
    bam_hdr_destroy(h); sam_close(in); bam_destroy1(b);
    return 0;
}

// trim polyA tail or polyT tail
int fxt_trim(int argc, char *argv[]) {
    if (argc != 6) {
        err_printf("Usage: fxtools trim in.fa/fq min_trim_length min_fraction window_size lead_nonA_len > out.fa\n");
        err_printf("\n");
        return 0;
    }
    gzFile infp = xzopen(argv[1], "r");
    size_t min_len = atoi(argv[2]), win_size = atoi(argv[4]), lead_nonA_len = atoi(argv[5]);
    float min_frac = atof(argv[3]);
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = stdout;

    int i, j, l, trim_s, trim_e;
    while (kseq_read(seq) >= 0)
    {
        trim_s = 0; trim_e = seq->seq.l-1;
        for (l = 0; l <= lead_nonA_len; ++l) {
            // polyT head
            for (i = l; i < seq->seq.l - win_size + 1; ++i) {
                int T_n = 0;
                for (j = i; j < i+win_size; ++j) {
                    if (seq->seq.s[j] == 'T' || seq->seq.s[j] == 'U' || seq->seq.s[j] == 't' || seq->seq.s[j] == 'u')
                        T_n++;
                }
                if (T_n / (win_size + 0.0) >= min_frac) {
                    // select the last 'T/U' in window
                    for (j = i+win_size-1; j >= i; --j) {
                        if (seq->seq.s[j] == 'T' || seq->seq.s[j] == 'U' || seq->seq.s[j] == 't' || seq->seq.s[j] == 'u') { 
                            trim_s = j+1;
                            break;
                        }
                    }
                } else break;
            }
            if (trim_s < min_len) {
                for (i = seq->seq.l - l - 1; i >= win_size - 1; --i) {
                    int A_n = 0;
                    for (j = i; j >= i-win_size+1; --j) {
                        if (seq->seq.s[j] == 'A' || seq->seq.s[j] == 'a')
                            A_n++;
                    }
                    if (A_n / (win_size + 0.0) >= min_frac) {
                        // select the first 'A' in window
                        for (j = i-win_size+1; j <= i; ++j) {
                            if (seq->seq.s[j]=='A' || seq->seq.s[j]=='a') {
                                trim_e = j-1;
                                break;
                            }
                        }
                    } else break;
                }
                if (seq->seq.l-1 - trim_e < min_len) trim_e = seq->seq.l-1;
            }
            if (trim_s > 0 || trim_e < seq->seq.l-1) break;
        }

        if (trim_s == 0 && trim_e == seq->seq.l-1) fprintf(outfp, ">%s", seq->name.s);
        else fprintf(outfp, ">%s_trim:%d_%d", seq->name.s, trim_s+1, trim_e+1);
        if (seq->comment.l > 0) fprintf(outfp, " %s", seq->comment.s);
        fprintf(outfp, "\n");

        for (i = trim_s; i <= trim_e; ++i) {
            fprintf(outfp, "%c", seq->seq.s[i]);
        }
        fprintf(outfp, "\n");
    }

    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

// trim polyA tail or polyT head
int fxt_trimF(int argc, char *argv[]) {
    if (argc != 6) {
        err_printf("Usage: fxtools trimF in.fa/fq min_trim_length min_fraction window_size lead_nonA_len > out.fa\n");
        err_printf("\n");
        return 0;
    }
    gzFile infp = xzopen(argv[1], "r");
    size_t min_len = atoi(argv[2]), win_size = atoi(argv[4]), lead_nonA_len = atoi(argv[5]);
    float min_frac = atof(argv[3]);
    kseq_t *seq;
    seq = kseq_init(infp);
    FILE *outfp = stdout;

    size_t i, j, l, trim_s, trim_e;
    while (kseq_read(seq) >= 0)
    {
        trim_s = 0; trim_e = seq->seq.l-1;
        // polyT head
        for (l = 0; l <= lead_nonA_len; ++l) {
            for (i = l; i < seq->seq.l - win_size + 1; ++i) {
                int T_n = 0;
                for (j = i; j < i+win_size; ++j) {
                    if (seq->seq.s[j] == 'T' || seq->seq.s[j] == 'U' || seq->seq.s[j] == 't' || seq->seq.s[j] == 'u')
                        T_n++;
                }
                if (T_n / (win_size + 0.0) >= min_frac) {
                    // select the last 'T/U' in window
                    for (j = i+win_size-1; j >= i; --j) {
                        if (seq->seq.s[j] == 'T' || seq->seq.s[j] == 'U' || seq->seq.s[j] == 't' || seq->seq.s[j] == 'u') {
                            trim_s = j+1;
                            break;
                        }
                    }
                } else break;
            }
            if (trim_s < min_len) {
                for (i = seq->seq.l - l - 1; i >= win_size - 1; --i) {
                    int A_n = 0;
                    for (j = i; j >= i-win_size+1; --j) {
                        if (seq->seq.s[j] == 'A' || seq->seq.s[j] == 'a')
                            A_n++;
                    }
                    if (A_n / (win_size + 0.0) >= min_frac) {
                        // select the first 'A' in window
                        for (j = i-win_size+1; j <= i; ++j) {
                            if (seq->seq.s[j]=='A' || seq->seq.s[j]=='a') {
                                trim_e = j-1;
                                break;
                            }
                        }
                    } else break;
                }
                if (seq->seq.l-1 - trim_e < min_len) trim_e = seq->seq.l-1;
            }
            if (trim_s > 0 || trim_e < seq->seq.l-1) break;
        }
        if (trim_s > 0 || trim_e < seq->seq.l-1) {
            fprintf(outfp, ">%s_trim:%ld_%ld", seq->name.s, trim_s+1, trim_e+1);
            if (seq->comment.l > 0) fprintf(outfp, " %s", seq->comment.s);
            fprintf(outfp, "\n");

            for (i = trim_s; i <= trim_e; ++i) {
                fprintf(outfp, "%c", seq->seq.s[i]);
            }
            fprintf(outfp, "\n");
        }
    }

    err_gzclose(infp);
    kseq_destroy(seq);
    err_fclose(outfp);
    return 0;
}

// add bam tag
int fxt_bam_add_seq(int argc, char *argv[]) {
    if (argc != 3) {
        err_printf("Usage: fxtools bam-add-seq in.bam tag > out.bam\n");
        err_printf("\n");
        return 0;
    }
    char *in_bam_fn = argv[1], *tag = argv[2];
    htsFile *out = hts_open("-", "wb");

    samFile *in; bam_hdr_t *hdr; bam1_t *rec = bam_init1();

    if ((in = sam_open(in_bam_fn, "rb")) == NULL) err_fatal(__func__, "fail to open \"%s\"\n", in_bam_fn);
    if ((hdr = sam_hdr_read(in)) == NULL) err_fatal(__func__, "fail to read header for \"%s\"\n", in_bam_fn);
    if (sam_hdr_write(out, hdr) < 0) err_fatal_simple("Fail to write bam header.");

    int r, k, seq_len; uint8_t *bam_bseq;
    while ((r = sam_read1(in, hdr, rec)) >= 0) {
        if (rec->core.flag & BAM_FUNMAP) continue;
        bam_bseq = bam_get_seq(rec);
        seq_len = rec->core.l_qseq;
        char *polya_seq = (char*)_err_malloc(sizeof(char) * (seq_len+1));
        if ((rec->core.flag & BAM_FREVERSE) != 0) { // reverse
            for (k = 0; k < seq_len; ++k) polya_seq[k] = bam_rc_qseq_table[bam_seqi(bam_bseq, seq_len-k-1)];
        } else { // forware
            for (k = 0; k < seq_len; ++k) polya_seq[k] = bam_qseq_table[bam_seqi(bam_bseq, k)];
        }
        polya_seq[k] = '\0';
        bam_aux_append(rec, tag, 'Z', seq_len+1, (uint8_t*)polya_seq);
        free(polya_seq);

        if (sam_write1(out, hdr, rec) < 0) err_fatal_simple("Fail to write bam record.");
    }

    if (hdr) bam_hdr_destroy(hdr); if (in) sam_close(in); if (rec) bam_destroy1(rec);
    hts_close(out);
    return 0;
}

int main(int argc, char*argv[]) {
    if (argc < 2) return usage();
    if (strcmp(argv[1], "filter") == 0 || strcmp(argv[1], "fl") == 0) fxt_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "filter-name") == 0 || strcmp(argv[1], "fn") == 0) fxt_filter_name(argc-1, argv+1);
    else if (strcmp(argv[1], "filter-qual") == 0 || strcmp(argv[1], "fq") == 0) fxt_filter_qual(argc-1, argv+1);
    else if (strcmp(argv[1], "filter-bam") == 0 || strcmp(argv[1], "fb") == 0) fxt_filter_bam(argc-1, argv+1);
    else if (strcmp(argv[1], "filter-bam-name") == 0 || strcmp(argv[1], "fbn") == 0) fxt_filter_bam_name(argc-1, argv+1);
    else if (strcmp(argv[1], "split-fx") == 0 || strcmp(argv[1], "sx") == 0) fxt_split_fx(argc-1, argv+1);
    else if (strcmp(argv[1], "fq2fa") == 0 || strcmp(argv[1], "qa") == 0) fxt_fq2fa(argc-1, argv+1);
    else if (strcmp(argv[1], "fa2fq") == 0 || strcmp(argv[1], "aq") == 0) fxt_fa2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "re-co") == 0 || strcmp(argv[1], "rc") == 0) fxt_re_co(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2bed") == 0 || strcmp(argv[1], "bb") == 0) fxt_bam2bed(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2fx") == 0 || strcmp(argv[1], "bf") == 0) fxt_bam2fx(argc-1, argv+1);
    else if (strcmp(argv[1], "seq-display") == 0 || strcmp(argv[1], "sd") == 0) fxt_seq_dis(argc-1, argv+1);
    else if (strcmp(argv[1], "cigar-parse") == 0 || strcmp(argv[1], "cp") == 0) fxt_cigar_parse(argc-1, argv+1);
    else if (strcmp(argv[1], "length-parse") == 0 || strcmp(argv[1], "lp") == 0) fxt_len_parse(argc-1, argv+1);
    else if (strcmp(argv[1], "merge-fa") == 0 || strcmp(argv[1], "mf") == 0) fxt_merge_fa(argc-1, argv+1);
    else if (strcmp(argv[1], "merge-filter-fa") == 0 || strcmp(argv[1], "mff") == 0) fxt_merge_filter_fa(argc-1, argv+1);
    else if (strcmp(argv[1], "duplicate-seq") == 0 || strcmp(argv[1], "ds") == 0) fxt_duplicate_seq(argc-1, argv+1);
    else if (strcmp(argv[1], "duplicate-read") == 0 || strcmp(argv[1], "dd") == 0) fxt_duplicate_read(argc-1, argv+1);
    else if (strcmp(argv[1], "error-parse") == 0 || strcmp(argv[1], "ep") == 0) fxt_error_parse(argc-1, argv+1);
    else if (strcmp(argv[1], "dna2rna") == 0 || strcmp(argv[1], "dr") == 0) fxt_dna2rna(argc-1, argv+1);
    else if (strcmp(argv[1], "rna2dna") == 0 || strcmp(argv[1], "rd") == 0) fxt_rna2dna(argc-1, argv+1);
    else if (strcmp(argv[1], "trim") == 0 || strcmp(argv[1], "tr") == 0) fxt_trim(argc-1, argv+1);
    else if (strcmp(argv[1], "trimF") == 0 || strcmp(argv[1], "tr") == 0) fxt_trimF(argc-1, argv+1);
    else if (strcmp(argv[1], "sam-flag") == 0 || strcmp(argv[1], "sf") == 0) fxt_sam_flag(argc-1, argv+1);
    else if (strcmp(argv[1], "bam-add-seq") == 0 || strcmp(argv[1], "bt") == 0) fxt_bam_add_seq(argc-1, argv+1);
    else {fprintf(stderr, "unknow command [%s].\n", argv[1]); return 1; }

    return 0;
}
