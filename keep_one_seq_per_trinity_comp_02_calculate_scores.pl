#!/usr/local/bin/perl -w

use warnings;
use strict;
#use List::MoreUtils qw(uniq);

if( scalar ( @ARGV ) != 2 ){
	my @temp=split(/\//, $0);
	die
"

USAGE: $temp[-1] <Raw data for each contig> <file containing the weights>

This script reads the output file of the previous script ('*print_data.pl')
and calculates the score for each Trinity contig, based on certain weights.

";
}

# load the weights
my %weights = ();

open ( WEIGHTS, $ARGV[1] ) or die;
while ( my $line = <WEIGHTS> ) {
	chomp ($line);
	next if $line =~ /^#/;
	next if $line =~ /^$/;

	my @f = split ( /\s*=\s*/, $line );
	$weights{$f[0]} = $f[1];
}
close (WEIGHTS);

# load the contig data
open ( CONTIG_DATA, $ARGV[0] ) or die;
while ( my $line = <CONTIG_DATA> ) {
	chomp ($line);
	my @ctg_data = split ( /\t/, $line );
	
	my %ctg_data = ();
	
	foreach my $one_feature (@ctg_data) {
		my @f = split ( /:/, $one_feature );
		$ctg_data{$f[0]} = $f[1];
	}
	
	# also add the score and score_string variables
	$ctg_data{score} = 0;
	$ctg_data{score_string} = "";
	
	# if the contig is a reverse complement of another contig, then its score
	# will be zero
	if ( $ctg_data{rc_contig} == 1 ) {
		$ctg_data{score_string} .= 'rc_contig';
		print $ctg_data{ctg_name}, "\t";
		print $ctg_data{length}, "\t";
		print $ctg_data{score}, "\t";
		print $ctg_data{score_string}, "\n";
		next; # don't calculate anything else for such contigs
	}
	
	# Start calculating the score
	# A) Based on the contig sequence
	
	# Against the Uniref100 database ###################################
	# reward from the proportion of the database sequence covered in the Uniref blastx
	my $one_score = $weights{db_proportion_cov_uniref} * $ctg_data{db_proportion_cov_uniref};
	$one_score = sprintf ( "%.1f", $one_score );
	$ctg_data{score} += $one_score;
	$ctg_data{score_string} .= "db_proportion_cov_uniref: ($weights{db_proportion_cov_uniref} * $ctg_data{db_proportion_cov_uniref}) = $one_score ## ";
	
	# reward from the proportion of the Trinity contig sequence covered in the Uniref blastx
	$one_score = $weights{ctg_proportion_cov_uniref} * $ctg_data{ctg_proportion_cov_uniref};
	$one_score = sprintf ( "%.1f", $one_score );
	$ctg_data{score} += $one_score;
	$ctg_data{score_string} .= "ctg_proportion_cov_uniref: ($weights{ctg_proportion_cov_uniref} * $ctg_data{ctg_proportion_cov_uniref}) = $one_score ## ";
	
	# reward from the length of the database sequcne covered in the Uniref blastx
	if ( $ctg_data{db_length_cov_uniref_max} != 0 ) {
		$one_score = $weights{db_length_cov_uniref} * $ctg_data{db_length_cov_uniref} / $ctg_data{db_length_cov_uniref_max};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "db_length_cov_uniref: $weights{db_length_cov_uniref} * $ctg_data{db_length_cov_uniref} / $ctg_data{db_length_cov_uniref_max} = $one_score ## ";
	}
	
	# reward from the length of the Trinity contig sequence covered in the Uniref blastx
	if ( $ctg_data{ctg_length_cov_uniref_max} != 0 ) {
		$one_score = $weights{ctg_length_cov_uniref} * $ctg_data{ctg_length_cov_uniref} / $ctg_data{ctg_length_cov_uniref_max};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "ctg_length_cov_uniref: $weights{ctg_length_cov_uniref} * $ctg_data{ctg_length_cov_uniref} / $ctg_data{ctg_length_cov_uniref_max} = $one_score ## ";
	}
	
	# reward/penalty from the strand of the Uniref match, relative to the read1/2 ratio
	if ( $ctg_data{log2ratio_max} != 0 ) {
		$one_score = $weights{strand_uniref} * $ctg_data{strand_uniref} * ( (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} );
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "strand_uniref: $weights{strand_uniref} * $ctg_data{strand_uniref} * ( (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} ) = $one_score ## ";
	}
	####################################################################


	# Against the nr database ##########################################
	# reward from the proportion of the database sequence covered in the nr blastx
	$one_score = $weights{db_proportion_cov_nr} * $ctg_data{db_proportion_cov_nr};
	$one_score = sprintf ( "%.1f", $one_score );
	$ctg_data{score} += $one_score;
	$ctg_data{score_string} .= "db_proportion_cov_nr: ($weights{db_proportion_cov_nr} * $ctg_data{db_proportion_cov_nr}) = $one_score ## ";
	
	# reward from the proportion of the Trinity contig sequence covered in the nr blastx
	$one_score = $weights{ctg_proportion_cov_nr} * $ctg_data{ctg_proportion_cov_nr};
	$one_score = sprintf ( "%.1f", $one_score );
	$ctg_data{score} += $one_score;
	$ctg_data{score_string} .= "ctg_proportion_cov_nr: ($weights{ctg_proportion_cov_nr} * $ctg_data{ctg_proportion_cov_nr}) = $one_score ## ";
	
	# reward from the length of the database sequcne covered in the nr blastx
	if ( $ctg_data{db_length_cov_nr_max} != 0 ) {
		$one_score = $weights{db_length_cov_nr} * $ctg_data{db_length_cov_nr} / $ctg_data{db_length_cov_nr_max};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "db_length_cov_nr: $weights{db_length_cov_nr} * $ctg_data{db_length_cov_nr} / $ctg_data{db_length_cov_nr_max} = $one_score ## ";
	}
	
	# reward from the length of the Trinity contig sequence covered in the nr blastx
	if ( $ctg_data{ctg_length_cov_nr_max} != 0 ) {
		$one_score = $weights{ctg_length_cov_nr} * $ctg_data{ctg_length_cov_nr} / $ctg_data{ctg_length_cov_nr_max};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "ctg_length_cov_nr: $weights{ctg_length_cov_nr} * $ctg_data{ctg_length_cov_nr} / $ctg_data{ctg_length_cov_nr_max} = $one_score ## ";
	}
	
	# reward/penalty from the strand of the nr match, relative to the read1/2 ratio
	if ( $ctg_data{log2ratio_max} != 0 ) {
		$one_score = $weights{strand_nr} * $ctg_data{strand_nr} * ( (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} );
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "strand_nr: $weights{strand_nr} * $ctg_data{strand_nr} * ( (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} ) = $one_score ## ";
	}
	####################################################################
	
	# Other scores based on the contig sequence ########################
	# penalty for having low sequencing coverage (<2x)
	if ( $ctg_data{seq_cov} < 2 ) {
		$one_score = $weights{seq_cov_2x} * ( $ctg_data{seq_cov} - 2 ) / 2;
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "seq_cov_2x: $weights{seq_cov_2x} * ( $ctg_data{seq_cov} - 2 ) / 2 = $one_score ## ";
	}
	
	# reward or penalty for having a negative or positive log2(read1/2 ratio), respectively
	if ( $ctg_data{log2ratio_max} != 0 ) {
		$one_score = $weights{log2ratio} * (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "log2ratio: $weights{log2ratio} * (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} = $one_score ## ";
	}
	
	# reward for the length of the contig ##############################
	$one_score = $weights{relative_length} * $ctg_data{length} / $ctg_data{length_max};
	$one_score = sprintf ( "%.1f", $one_score );
	$ctg_data{score} += $one_score;
	$ctg_data{score_string} .= "relative_length: $weights{relative_length} * $ctg_data{length} / $ctg_data{length_max} = $one_score ## ";
	####################################################################
	
	
	# B) Based on the called ORFs ######################################
	# penalty if >1 ORFs are called
	if ( $ctg_data{orf_number} > 0 ) {
		$one_score = $weights{orf_number} * ( 1 - $ctg_data{orf_number} ) / $ctg_data{orf_number};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "orf_number: $weights{orf_number} * ( 1 - $ctg_data{orf_number} ) / $ctg_data{orf_number} = $one_score ## ";
	}
	
	# penalty if each ORF matches a different protein in nr
	if ( $ctg_data{number_of_different_matches} > 0 ) {
		$one_score = $weights{number_of_different_matches} * ( 1 - $ctg_data{number_of_different_matches} ) / $ctg_data{number_of_different_matches};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "number_of_different_matches: $weights{number_of_different_matches} * ( 1 - $ctg_data{number_of_different_matches} ) / $ctg_data{number_of_different_matches} = $one_score ## ";
	}
	
	# reward or penalty depending on which strand the ORFs are called (in
	# the case where all ORFs are called in one strand)
	if ( $ctg_data{log2ratio_max} != 0 ) {
		if ( $ctg_data{number_of_plus_ORFs} > 0 && $ctg_data{number_of_minus_ORFs} == 0 ) {
			my $orf_strand = 1;
			$one_score = $weights{orfs_strand} * $orf_strand * (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max};
			$one_score = sprintf ( "%.1f", $one_score );
			$ctg_data{score} += $one_score;
			$ctg_data{score_string} .= "orf_strand: $weights{orfs_strand} * $orf_strand * (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} = $one_score ## ";
		}
		elsif ( $ctg_data{number_of_plus_ORFs} == 0 && $ctg_data{number_of_minus_ORFs} > 0 ) {
			my $orf_strand = -1;
			$one_score = $weights{orfs_strand} * $orf_strand * (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max};
			$one_score = sprintf ( "%.1f", $one_score );
			$ctg_data{score} += $one_score;
			$ctg_data{score_string} .= "orf_strand: $weights{orfs_strand} * $orf_strand * (-1) * $ctg_data{log2ratio} / $ctg_data{log2ratio_max} = $one_score ## ";
		}
	}
	
	# penalty if ORFs are called in both strands. Moreover, if any of them,
	# or both of them have a nr match, the penalty increases
	if ( $ctg_data{number_of_plus_ORFs} > 0 && $ctg_data{number_of_minus_ORFs} > 0 ) {
		# pick the right weight depending on whether there are matches
		my $weight_tmp;
		if ( $ctg_data{number_of_plus_ORFs_with_match} > 0 && $ctg_data{number_of_minus_ORFs_with_match} > 0 ) {
			# ORFs in both strands have match(es).
			$weight_tmp = $weights{orfs_both_strands_both_with_match};
		}
		elsif ( $ctg_data{number_of_plus_ORFs_with_match} > 0 && $ctg_data{number_of_minus_ORFs_with_match} == 0 ) {
			$weight_tmp = $weights{orfs_both_strands_one_with_match};
		}
		elsif ( $ctg_data{number_of_plus_ORFs_with_match} == 0 && $ctg_data{number_of_minus_ORFs_with_match} > 0 ) {
			$weight_tmp = $weights{orfs_both_strands_one_with_match};
		}
		elsif ( $ctg_data{number_of_plus_ORFs_with_match} == 0 && $ctg_data{number_of_minus_ORFs_with_match} == 0 ) {
			$weight_tmp = $weights{orfs_both_strands_none_with_match};
		}
		#else {
		#	print STDERR "You're not supposed to see this!\n";
		#}
		
		# now calculate the score
		my @tmp = ();
		push ( @tmp, $ctg_data{number_of_plus_ORFs} );
		push ( @tmp, $ctg_data{number_of_minus_ORFs} );
		
		@tmp = sort { $a <=> $b } @tmp;
		
		$one_score = $weight_tmp * (-1) * $tmp[0] / $tmp[1];
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "orfs_both_strands: $weight_tmp * (-1) * $tmp[0] / $tmp[1] = $one_score ## ";
	}
	####################################################################
	
	# C) Based on the sequencing coverage dips #########################
	# penalty depending on the number of sequencing coverage dips
	if ( $ctg_data{seq_cov_dips_max} != 0 ) {
		$one_score = $weights{seq_cov_dips} * (-1) * $ctg_data{seq_cov_dips} / $ctg_data{seq_cov_dips_max};
		$one_score = sprintf ( "%.1f", $one_score );
		$ctg_data{score} += $one_score;
		$ctg_data{score_string} .= "seq_cov_dips: $weights{seq_cov_dips} * (-1) * $ctg_data{seq_cov_dips} / $ctg_data{seq_cov_dips_max} = $one_score ## ";
	}
	####################################################################
	
	
	# print the contig name, score and score_string
	print $ctg_data{ctg_name}, "\t";
	print $ctg_data{length}, "\t";
	print sprintf ( "%.1f", $ctg_data{score} ), "\t";
	print $ctg_data{score_string}, "\n";
}
