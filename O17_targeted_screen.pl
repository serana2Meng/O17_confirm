#!/usr/bin/env perl
#
# O17_targeted_screen.pl - O17:K4 serotype confirmation script
#
# Purpose: Confirm O17 candidates from BLAST m6 results.
#          Assumes BLAST was run with genome as query, reference sequences as subject.
#          Two-tier strategy: ① Full module BLAST (coverage ≥90% → O17)
#                             ② Individual gene search + collinearity check (fallback)
#
# Author: [Xiaohua Meng]
#
# ============================================================

use strict;
use warnings;
use Getopt::Long;

# ---------- Help message ----------
sub print_help {
    print <<'EOF';

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

THRESHOLDS (currently set in script):
  Full module coverage threshold:  90%
  Single gene coverage threshold:  90%
  Maximum collinearity gap:        5000 bp

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

EOF
    exit 0;
}

# ---------- Parse command-line arguments ----------
my ($sample_id, $infile, $outfile, $help);
GetOptions(
    "s=s" => \$sample_id,
    "b=s" => \$infile,
    "o=s" => \$outfile,
    "h|help" => \$help,
) or die "Error: Invalid options. Use -h for help.\n";

if ($help || !($sample_id && $infile && $outfile)) {
    print_help();
}

die "Error: Sample ID (-s) is required!\n" unless $sample_id;
die "Error: BLAST file '$infile' not found!\n" unless -f $infile;
die "Error: Output file (-o) is required!\n" unless $outfile;

# ---------- Thresholds ----------
my $MODULE_COV_THRESHOLD = 0.90;   # Minimum coverage for full‑module hit
my $GENE_COV_THRESHOLD   = 0.90;   # Minimum coverage for individual genes
my $COLLINEARITY_MAX_GAP = 5000;   # Maximum allowed gap between collinear genes (bp)

# ---------- Reference sequence lengths (must match your actual data) ----------
my %REF_LEN = (
    'module_full' => 3531,
    'gene_273'    => 921,
    'gene_274'    => 804,
    'gene_275'    => 1668,
);

# ---------- Read and parse m6 file ----------
# m6 format: qseqid sseqid pident length qstart qend sstart send evalue bitscore
# Here: qseqid = contig, sseqid = reference sequence ID
my %hits;  # $hits{$base} = [ { hit1 }, { hit2 }, ... ]
open(BLAST, "<", $infile) or die "Cannot open $infile: $!";
while (<BLAST>) {
    chomp;
    my ($qseqid, $sseqid, $pident, $length, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split /\t/;
    # Extract base name from subject ID (sseqid) because reference sequences are in subject
    my ($base) = $sseqid =~ /^(module_full|gene_273|gene_274|gene_275)/;
    next unless $base;
    push @{ $hits{$base} }, {
        contig   => $qseqid,   # contig name (query)
        start    => $qstart < $qend ? $qstart : $qend,  # genomic start (query)
        end      => $qstart < $qend ? $qend : $qstart,  # genomic end (query)
        strand   => $qstart < $qend ? '+' : '-',        # orientation from query coordinates
        qseqid   => $qseqid,
        sseqid   => $sseqid,
        pident   => $pident,
        length   => $length,
        qstart   => $qstart,
        qend     => $qend,
        sstart   => $sstart,
        send     => $send,
        evalue   => $evalue,
        bitscore => $bitscore,
    };
}
close BLAST;

# ---------- Extract the best hit for each query (by coverage) ----------
sub best_hit {
    my ($base) = @_;
    my $hits_ref = $hits{$base};
    return undef unless $hits_ref && @$hits_ref;
    my $best = undef;
    my $best_cov = 0;
    my $qlen = $REF_LEN{$base};
    foreach my $hit (@$hits_ref) {
        my $hit_len = $hit->{length};   # BLAST alignment length
        my $cov = $hit_len / $qlen;
        if ($cov > $best_cov) {
            $best_cov = $cov;
            $best = $hit;
            $best->{coverage} = $cov;
        }
    }
    return $best;
}

my $module_best = best_hit('module_full');
my $best_273    = best_hit('gene_273');
my $best_274    = best_hit('gene_274');
my $best_275    = best_hit('gene_275');

# ---------- Helper: extract genomic coordinates from a hit ----------
sub get_coords {
    my ($hit) = @_;
    return undef unless $hit;
    return {
        contig   => $hit->{contig},
        start    => $hit->{start},
        end      => $hit->{end},
        strand   => $hit->{strand},
        coverage => $hit->{coverage},
        pident   => $hit->{pident},
        length   => $hit->{length},
    };
}

my $coords_273 = get_coords($best_273);
my $coords_274 = get_coords($best_274);
my $coords_275 = get_coords($best_275);

# ---------- Decision variables ----------
my $result_code = '';      # 'O17' or 'O17*'
my $should_print = 0;

# ---------- ROUND 1: Full‑module confirmation (primary criterion) ----------
if ($module_best && $module_best->{coverage} >= $MODULE_COV_THRESHOLD) {
    $result_code = 'O17';
    $should_print = 1;
    goto OUTPUT;
}

# ---------- ROUND 2: Individual gene search + collinearity check (fallback) ----------
my $all_present = 1;
my $missing = [];

foreach my $g (['gene_273', $coords_273], ['gene_274', $coords_274], ['gene_275', $coords_275]) {
    my ($name, $c) = @$g;
    unless ($c && $c->{coverage} >= $GENE_COV_THRESHOLD) {
        $all_present = 0;
        push @$missing, $name;
    }
}

if (!$all_present) {
    $should_print = 0;
    goto OUTPUT;
}

# Check if all three genes are on the same contig
my $contig_273 = $coords_273->{contig};
my $contig_274 = $coords_274->{contig};
my $contig_275 = $coords_275->{contig};

if ($contig_273 eq $contig_274 && $contig_274 eq $contig_275) {
    # Same contig: verify order and gaps
    my @genes = ($coords_273, $coords_274, $coords_275);
    @genes = sort { $a->{start} <=> $b->{start} } @genes;
    my @ids = map {
        $_ == $coords_273 ? '273' :
        ($_ == $coords_274 ? '274' : '275')
    } @genes;
    my $order_ok = ($ids[0] eq '273' && $ids[1] eq '274' && $ids[2] eq '275');
    my $gap_ok = 1;
    for (my $i = 0; $i < @genes - 1; $i++) {
        my $gap = $genes[$i+1]->{start} - $genes[$i]->{end};
        if ($gap > $COLLINEARITY_MAX_GAP) { $gap_ok = 0; last; }
    }
    if ($order_ok && $gap_ok) {
        $result_code = 'O17';
        $should_print = 1;
    } else {
        $should_print = 0;
    }
} else {
    # Genes on different contigs → candidate (O17*)
    $result_code = 'O17*';
    $should_print = 1;
}

# ---------- Output ----------
OUTPUT:
if ($should_print) {
    open(OUT, ">", $outfile) or die "Cannot write to $outfile: $!";
    print OUT $sample_id . "\t" . $result_code . "\n";
    close OUT;
} else {
    open(OUT, ">", $outfile) or die "Cannot write to $outfile: $!";
    close OUT;  # empty file
}
