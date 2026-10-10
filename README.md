# O17_targeted_screen.pl

A Perl script for **downstream confirmation of O17:K4 serotype candidates** from BLASTn m6 results.

This script is designed to run after an initial O17 screening step.  
It takes samples that were either:

- flagged as **O17-positive** by the initial screen, or
- flagged as **requiring manual review** (e.g., include `O17` candidates),

and performs a targeted check to confirm whether the O17-specific region is intact or disrupted.

It uses a two-tier strategy:

1. **Full module BLAST** — if the full O17 module is covered ≥ 90%, the sample is confirmed as O17.
2. **Individual gene search and collinearity check** (fallback) — if the full module is not covered, the script checks three specific genes (`gene_273`, `gene_274`, `gene_275`). If all three are present at ≥ 90% coverage and arranged collinearly on the same contig, the sample is confirmed as O17. If all three genes are present but on different contigs, the sample is reported as an O17 candidate (`O17*`) for manual review.

## Usage

```bash
perl O17_targeted_screen.pl -s <sample_id> -b <blast.m6> -o <output.txt>
```

### Full help message

```
O17_targeted_screen.pl - O17:K4 serotype confirmation script

DESCRIPTION:
  This script confirms O17 candidates from BLAST m6 output.
  It uses a two-tier strategy:
    1) Full module BLAST (coverage ≥90% → O17 positive)
    2) Individual gene search + collinearity check (fallback mode)

USAGE:
  perl O17_targeted_screen.pl -s <sample_id> -b <blast.m6> -o <output.txt>

REQUIRED ARGUMENTS:
  -s <sample_id>   Sample identifier (used in output, e.g., GCA_123456)
  -b <blast.m6>    BLAST m6 output file (from blastn -outfmt 6)
  -o <output.txt>  Output file (one line if positive/candidate, empty if negative)

OPTIONAL:
  -h, --help       Show this help message and exit

OUTPUT FORMAT:
  If O17 positive (confirmed by full module or collinearity):
    <sample_id>\tO17
  If O17 candidate (genes exist but on different contigs):
    <sample_id>\tO17*
  If negative:
    (empty file, no output)

NOTES:
  - BLAST must be run with genome as query and reference sequences as subject.
  - The script expects query sequence IDs to start with these names
    (suffix like ':length' is allowed, e.g., gene_273:921).
```
