fxtools
=======
 ```
Program: fxtools (light-weight processing tool for FASTA, FASTQ and BAM format data)
Usage:   fxtools <command> [options]

Command: 
         filter (fl)           filter fa/fq sequences with specified length boundary.
         filter-name (fn)      filter fa/fq sequences with specified name.
         filter-bam (fb)       filter bam/sam records with specified read length boundary.
         filter-bam-name (fbn) filter bam/sam records with specified read name.
         fq2fa (qa)            convert FASTQ format data to FASTA format data.
         fa2fq (aq)            convert FASTA format data to FASTQ format data.
         bam2bed (bb)          convert BAM file to BED file. seperated exon regions for spliced BAM
         re-co (rc)            convert DNA sequence(fa/fq) to its reverse-complementary sequence.
         seq-display (sd)      display a specified region of FASTA/FASTQ file.
         cigar-parse (cp)      parse the given cigar(stdout).
         length-parse (lp)     parse the length of sequences in fa/fq file.
         merge-fa (mf)         merge the reads with same read name in fasta/fastq file.
         merge-filter-fa (mff) merge and filter the reads with same read name in fasta file.
         error-parse (ep)      parse indel and mismatch error based on CIGAR and NM in bam file.
         dna2rna (dr)          convert DNA fa/fq to RNA fa/fq.
         rna2dna (rd)          convert RNA fa/fq to DNA fa/fq.
```
