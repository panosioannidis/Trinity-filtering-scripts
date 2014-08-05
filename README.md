Trinity-filtering-scripts
=========================

The scripts for filtering of the raw Trinity output transcript/contig set

I will keep adding information on what each script does and how to link them together

1\. <b>keep_one_seq_per_trinity_comp_01_print_data.pl</b>

{{{USAGE: keep_one_seq_per_trinity_comp_01_print_data.pl \<file containing paths to input files\>}}}

Prints data for each contig, such as its length, features of the similarities to other sequences, number of sequencing coverage dips. It outputs a tab-delimited file and the data for each contig are printed in one line, in a "key:value" format.

2\. <b>keep_one_seq_per_trinity_comp_02_calculate_scores.pl</b>

USAGE: keep_one_seq_per_trinity_comp_02_calculate_scores.pl \<Raw data for each contig\> \<file containing the weights\>

This script calculates the score for each contig using the output of the previous script and also by taking into account the weights for certain features.

3\. <b>keep_one_seq_per_trinity_comp_03_select_contigs.pl</b>

USAGE: keep_one_seq_per_trinity_comp_03_select_contigs.pl \<List of contigs with their scores\> \<score cutoff\>

It selects up to one contig per Trinity component, based on the score calculated in the previous step.
