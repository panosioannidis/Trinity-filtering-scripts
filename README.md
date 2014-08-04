Trinity-filtering-scripts
=========================

The scripts for filtering of the raw Trinity output transcript/contig set

I will keep adding information on what each script does and how to link them together

1. 01_print_data.pl

Prints data for each contig, such as its length, features of the similarities to other sequences, number of sequencing coverage dips. It outputs a tab-delimited file and the data are printed in a "key:value" format.

2. 02_calculate_scores.pl

This script calculates the score for each contig by taking into account the weights for certain features (present in a separate file).

3. 03_select_contigs.pl

It selects up to one contig per Trinity component, based on the score.
