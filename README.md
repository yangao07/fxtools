# fxtools
light-weight processing tool for FASTA, FASTQ and BAM format data

## Installation
```
git clone https://github.com/yangao07/fxtools.git --recursive
cd fxtools; make
```

## Commands and options
```
Program: fxtools (light-weight processing tool for FASTA, FASTQ and BAM format data)
Usage:   fxtools <command> [options]

Command:
         filter (fl)           filter fa/fq sequences with specified length boundary.
         filter-name (fn)      filter fa/fq sequences with specified name.
         filter-qual (fq)      filter fq sequences with specified min. average quality.
         filter-bam (fb)       filter bam/sam records with specified read length boundary.
         filter-bam-name (fbn) filter bam/sam records with specified read name.
         split-fx (sx)         split fa/fq file into multipule files.
         fq2fa (qa)            convert FASTQ format data to FASTA format data.
         fa2fq (aq)            convert FASTA format data to FASTQ format data.
         bam2bed (bb)          convert BAM file to BED file. seperated exon regions for spliced BAM
         bam2fx (bf)           convert BAM file to FASTA/FASTQ file.
         re-co (rc)            convert DNA sequence(fa/fq) to its reverse-complementary sequence.
         seq-display (sd)      display a specified region of FASTA/FASTQ file.
         cigar-parse (cp)      parse the given cigar(stdout).
         length-parse (lp)     parse the length of sequences in fa/fq file.
         merge-fa (mf)         merge the reads with same read name in a fasta/fastq file.
         merge-fas (mfs)       merge the reads with same read name in two fasta/fastq files.
         merge-filter-fa (mff) merge and filter the reads with same read name in fasta file.
         duplicate-seq (ds)    duplicate all the read sequences with specific copy number.
         duplicate-read (dd)   duplicate all the read records with specific copy number.
         error-parse (ep)      parse indel and mismatch error based on CIGAR and NM in SAM/BAM/GAF file.
         dna2rna (dr)          convert DNA fa/fq to RNA fa/fq.
         rna2dna (rd)          convert RNA fa/fq to DNA fa/fq.
         trim (tr)             trim polyA tail(polyT head).
         trimF (tf)            trim and filter with polyA tail(polyT head). Only polyA reads will be kept.
         bam-add-seq (bt)      add polyA sequence to BAM record.
```
