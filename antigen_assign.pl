#!/usr/bin/perl
use strict;
use warnings;

die "Usage:\nperl $0 <in> <out> <percent> <ratio> <gene> <name>\n" unless @ARGV==6;
my($input,$out,$per,$ra,$list,$name)=@ARGV[0,1,2,3,4,5];

my %hash;
my %fragments;
my %gene_len;
my %contig_len;
my $auto_done = 0;
my $decision_type = "auto";
my @candidates_info;

# ==================== 读取contig长度 ====================
open IN,"<$input" or die $!;
while(<IN>){
    chomp;
    my @ar=split/\t/,$_;
    my $contig = $ar[0];
    if ($contig =~ /length_(\d+)/) {
        $contig_len{$contig} = $1;
    }
}
close IN;

# 读取基因列表
open NAME,"<$list" or die $!;
while(<NAME>){
    chomp;
    $hash{$_}=0;
}
close NAME;

# ==================== 第1级：找完整匹配 ====================
open IN,"<$input" or die $!;
while(<IN>){
    chomp;
    my @ar=split/\t/,$_;
    my @s=split/\:/,$ar[1];
    my $gene = $s[0];
    my $gene_len = $s[1];
    $gene_len{$gene} = $gene_len;
    my $ali_len = abs($ar[9]-$ar[8])+1;
    my $ratio = $ali_len / $gene_len;
    
    if($ar[2] >= $per && $ratio >= $ra){
        $hash{$gene}=1;
        $auto_done = 1;
    }
}
close IN;

# ==================== 片段收集 ====================
open IN,"<$input" or die $!;
while(<IN>){
    chomp;
    my @ar=split/\t/,$_;
    my @s=split/\:/,$ar[1];
    my $gene = $s[0];
    my $gene_len = $s[1];
    $gene_len{$gene} = $gene_len;
    my $ali_len = abs($ar[9]-$ar[8])+1;
    
    if($ar[2] >= 80 && $ali_len > 200){
        my $ref_start = $ar[8] < $ar[9] ? $ar[8] : $ar[9];
        my $ref_end   = $ar[8] > $ar[9] ? $ar[8] : $ar[9];
        my $q_start   = $ar[6] < $ar[7] ? $ar[6] : $ar[7];
        my $q_end     = $ar[6] > $ar[7] ? $ar[6] : $ar[7];
        my $q_id      = $ar[0];
        push @{$fragments{$gene}}, [$ref_start, $ref_end, $ali_len, $q_start, $q_end, $q_id];
    }
}
close IN;

# ==================== 屏幕打印函数 ====================
sub print_merged_stats {
    print "\n" . "=" x 80 . "\n";
    print "O-antigen merging statistics for $name:\n";
    print "=" x 80 . "\n";
    
    foreach my $gene (sort keys %fragments) {
        my @sorted = sort { $a->[0] <=> $b->[0] } @{$fragments{$gene}};
        
        my $total_len = 0;
        my $max_frag_len = 0;
        my ($cur_start, $cur_end) = ($sorted[0]->[0], $sorted[0]->[1]);
        $max_frag_len = $sorted[0]->[2] if $sorted[0]->[2] > $max_frag_len;
        
        my @merged_intervals;
        push @merged_intervals, "$cur_start-$cur_end";
        
        for (my $i = 1; $i < @sorted; $i++) {
            my ($s, $e, $len) = @{$sorted[$i]};
            $max_frag_len = $len if $len > $max_frag_len;
            
            if ($s <= $cur_end + 1) {
                if ($e > $cur_end) {
                    $cur_end = $e;
                    $merged_intervals[-1] = "$cur_start-$cur_end";
                }
            } else {
                $total_len += ($cur_end - $cur_start + 1);
                ($cur_start, $cur_end) = ($s, $e);
                push @merged_intervals, "$cur_start-$cur_end";
            }
        }
        $total_len += ($cur_end - $cur_start + 1);
        
        my $gene_len = $gene_len{$gene};
        my $cover_ratio = $total_len / $gene_len;
        my $max_frag_ratio = $max_frag_len / $gene_len;
        
        print "Gene: $gene (length: $gene_len)\n";
        print "  Merged intervals: " . join(", ", @merged_intervals) . "\n";
        print "  Total covered: $total_len bp (Cover ratio: " . sprintf("%.4f", $cover_ratio) . ")\n";
        print "  Max fragment: $max_frag_len bp (Ratio: " . sprintf("%.4f", $max_frag_ratio) . ")\n";
        print "\n";
    }
    print "=" x 80 . "\n\n";
}

# ==================== 完整匹配分支 ====================
if ($auto_done) {
    open OT,">$out" or die $!;
    foreach my $i (sort keys %hash){
        if ($hash{$i} == 1) {
            print OT "$name\t$i\t1\n";
        }
    }
    close OT;
    print_merged_stats();
    exit 0;
}

# ==================== 无完整匹配：合并区间，计算覆盖度 ====================
my @candidates;
foreach my $gene (keys %fragments){
    # 按contig分组合并（同一contig内去重，不同contig累加）
    my %by_contig;
    foreach my $f (@{$fragments{$gene}}) {
        my $q_id = $f->[5];
        push @{$by_contig{$q_id}}, $f;
    }
    
    my $total_len = 0;
    my $max_frag_len = 0;
    
    foreach my $contig (keys %by_contig) {
        my @sorted = sort { $a->[0] <=> $b->[0] } @{$by_contig{$contig}};
        
        my $contig_len = 0;
        my ($cur_start, $cur_end) = ($sorted[0]->[0], $sorted[0]->[1]);
        $max_frag_len = $sorted[0]->[2] if $sorted[0]->[2] > $max_frag_len;
        
        for (my $i = 1; $i < @sorted; $i++) {
            my ($s, $e, $len) = @{$sorted[$i]};
            $max_frag_len = $len if $len > $max_frag_len;
            
            if ($s <= $cur_end + 1) {
                $cur_end = $e if $e > $cur_end;
            } else {
                $contig_len += ($cur_end - $cur_start + 1);
                ($cur_start, $cur_end) = ($s, $e);
            }
        }
        $contig_len += ($cur_end - $cur_start + 1);
        $total_len += $contig_len;
    }
    
    my $gene_len = $gene_len{$gene};
    my $cover_ratio = $total_len / $gene_len;
    my $max_frag_ratio = $max_frag_len / $gene_len;
    
    if ($cover_ratio > 0.5) {
        push @candidates, [$gene, $cover_ratio, $max_frag_len, $max_frag_ratio];
    }
}

# ==================== 输出主文件 ====================
open OT,">$out" or die $!;

if (!@candidates) {
    close OT;
    $decision_type = "manual";
    @candidates_info = ();
    generate_list_file("manual");
    print_merged_stats();
    exit 0;
}

my @qualified = grep { $_->[1] > 0.8 } @candidates;

if (!@qualified) {
    close OT;
    $decision_type = "manual";
    @candidates_info = sort { $b->[1] <=> $a->[1] } @candidates;
    generate_list_file("manual");
    print_merged_stats();
    exit 0;
}

# ==================== 判定逻辑（无质量评分） ====================
# 首先按最长片段比例降序排序
@qualified = sort { $b->[3] <=> $a->[3] } @qualified;

my $can_auto = 1;
my $winner = undef;

if (@qualified == 1) {
    $winner = $qualified[0]->[0];
} else {
    my ($gene1, $ratio1, $max_len1, $max_ratio1) = @{$qualified[0]};
    my ($gene2, $ratio2, $max_len2, $max_ratio2) = @{$qualified[1]};
    
    # 1. 比较最长片段比例
    if ($max_ratio1 / $max_ratio2 > 1.5) {
        $winner = $gene1;
    } else {
        # 2. 最长片段比例接近，比较综合覆盖度
        if ($ratio1 - $ratio2 > 0.05) {
            $winner = $gene1;
        } else {
            # 3. 两者都接近，无法自动判定
            $can_auto = 0;
        }
    }
}

if (!$can_auto) {
    close OT;
    $decision_type = "manual";
    @candidates_info = @qualified;
    generate_list_file("manual");
    print_merged_stats();
    exit 0;
}

$hash{$winner} = "1*";
$decision_type = "merge";

foreach my $i (sort keys %hash){
    if ($hash{$i} eq "1*" || $hash{$i} == 1) {
        print OT "$name\t$i\t$hash{$i}\n";
    }
}
close OT;

@candidates_info = @qualified;
generate_list_file("merge");
print_merged_stats();

# ==================== 子函数 ====================
sub generate_list_file {
    my ($type) = @_;
    
    my $list_file = ($type eq "merge") ? "$out.merge" : "$out.manual";
    
    open LIST,">$list_file" or die $!;
    print LIST "Sample: $name\n";
    print LIST "Decision_Type: $type\n";
    print LIST "Input_File: $input\n";
    print LIST "-" x 70 . "\n";
    print LIST "Rank\tO_Antigen\tCover_Ratio\tMax_Fragment\tMax_Frag_Ratio\n";
    print LIST "-" x 70 . "\n";
    
    if (@candidates_info) {
        my @sorted = sort { $b->[1] <=> $a->[1] } @candidates_info;
        my $max_show = @sorted < 5 ? @sorted : 5;
        for (my $i = 0; $i < $max_show; $i++) {
            my ($gene, $ratio, $max_len, $max_ratio) = @{$sorted[$i]};
            printf LIST "%d\t%s\t%.4f\t%d\t%.4f\n", $i+1, $gene, $ratio, $max_len, $max_ratio;
        }
    } else {
        print LIST "No candidates found.\n";
    }
    
    print LIST "-" x 70 . "\n";
    if ($type eq "merge") {
        print LIST "Action: Please manually verify the top candidate(s).\n";
    } else {
        print LIST "Action: Cannot determine automatically. Please manually check.\n";
    }
    close LIST;
}
