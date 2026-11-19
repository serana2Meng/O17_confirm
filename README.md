# O17_targeted_screen.pl

A Perl script for **downstream confirmation of O17 serotype candidates** from BLASTn m6 results.

This script is designed to run after an initial O17 screening step. 
It takes samples that fall into one of the following three categories:

- flagged as **sample O17**   by the initial screen
- flagged as **sample O17***  O17 candidates
Only samples flagged as O17 or O17* by the initial screen are carried forward. Samples with no detectable O17 signal are not pursued further and produce no output.

and performs a targeted check to confirm whether the O17-specific region is intact or disrupted.

It uses a two-tier strategy:

- **Full module BLAST** — if the full O17 module is covered ≥ 90%, the sample is confirmed as O17.
- **Individual gene search and collinearity check** (fallback) — if the full module is not covered, the script checks three specific genes (`gene_273`, `gene_274`, `gene_275`). If all three are present at ≥ 90% coverage and arranged collinearly on the same contig, the sample is confirmed as O17. If all three genes are present but on different contigs, the sample is reported as an O17 candidate (`O17*`) for manual review.

## Dependencies

- Perl (>= 5.10)
- BLAST+ (makeblastdb, blastn)

## Example workflow

```bash
# Build BLAST database from the O17 target reference
makeblastdb -in O17_references_singleline.fasta -dbtype nucl -out O17_target_db

# Run BLASTn against the O17 target region

blastn \
  -db O17_target_db \
  -query candidate_genome.fna \
  -outfmt 6 \
  -evalue 1e-5 \
  -num_threads 4 \
  -out candidate_genome.blast.m6

# Run the targeted screening script on the BLAST output

perl O17_targeted_screen.pl -s candidate_genome -b candidate_genome.blast.m6 -o candidate_genome_O17.txt

The -evalue 1e-5 threshold is set at the BLAST step; the Perl script does not re-apply this threshold.

text
```
## Usage

```bash
perl O17_targeted_screen.pl -s <sample_id> -b <sample.blast.m6> -o <output.txt>
```
## Parameters
  -s <sample_id>   Sample identifier (used in output, e.g., GCA_123456)
  -b <blast.m6>    BLAST m6 output file (from blastn -outfmt 6)
  -o <output.txt>  Output file (one line if positive/candidate, empty if negative)

OPTIONAL:
  -h, --help       Show this help message and exit
    
See example/README.md for details.
