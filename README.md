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

## Dependencies
