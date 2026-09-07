## Script function
perl O17_targeted_screen.pl -s GCA_015798135.1_PDT000396168.1_genomic -b GCA_002GCA_015798135.1_PDT000396168.1_genomic.blast.m6 -o GCA_015798135.1_PDT000396168.1_genomic.result.txt
This script parses BLAST alignment results of target sequences to judge whether target sequences are present. It outputs the target‑sequence name when a target is confirmed present; outputs the target name followed by an asterisk (`*`) for probable presence; and produces no output for other cases.

perl /Genome/mengxh/Pathogen/Phylogenetic_trees/VP/Ref_O/antigen_assign_merge2.pl GCA_015798135.1.blast.m6 GCA_015798135.1.result 90 0.9 O1_O18.list GCA_015798135.1
This script is used for O-antigen serotype determination. It outputs the definite O-antigen serotype when the criteria are fully satisfied. If the sample indirectly meets the screening conditions, it returns the serotype with an asterisk suffix (serotype*). For undetermined cases, the script outputs basic statistical evidence for manual curation and further judgment.
