#!/usr/local/bin/perl -w

use warnings;
use strict;
use List::MoreUtils qw(uniq);

if( scalar ( @ARGV ) != 1 ){
	my @temp=split(/\//, $0);
	die
"

USAGE: $temp[-1] <file containing paths to input files>

This script loads evidence for each Trinity contig and prints it out in one
line per contig. The next script ('*_calculate_scores.pl') will take this
file and calculate the scores for each contig.

";
}

# open the file with the paths to input files
my %input = ();

open ( INPUT, $ARGV[0] ) or die;
while ( my $line = <INPUT> ) {
	chomp ($line);
	my @f = split ( /\s*=\s*/, $line );
	$input{$f[0]} = $f[1];
}
close (INPUT);

# Part A: Load the data ################################################

# A1) Load the contig and component names
my %contig = ();
my %comp_max = (); # will hold the maximum of certain values per component
LoadContigs ( $input{'trinity.fasta'} );

# A2) Remove contigs that are the reverse complement of other contigs
TagRCContigs ( $input{'rc_contigs'} );

# A3) Load the ORFs called for each contig
LoadORFs ( $input{'orf_calls'} );

# A4) add the log2(read12 ratio) for each contig
LoadLogRatiosSeqCov ($input{'log2_ratios'});

# A5) read the BLASTX of contigs against Uniref100 (*qr_sb_cov file)
LoadBlastxResults ( $input{'ctg_vs_uniref'}, "uniref" );

# A6) read the BLASTX of contigs against nr
LoadBlastxResults ( $input{'ctg_vs_nr'}, "nr" );

# A7) Load the BLASTP results of TransDecoder ORFs against nr
LoadOrfBlastpResults ( $input{'tdec_orfs_vs_nr'}, "tdec" );

# A8) Load the BLASTP results of GetOrf ORFs against nr
LoadOrfBlastpResults ( $input{'getorf_orfs_vs_nr'}, "getorf" );

# A9) Load the sequencing coverage dips for each contig
LoadSeqCovDips ( $input{'seq_cov_dips'} );

########################################################################


# Part B: Print the results out ########################################

my @ctg_names = keys(%contig);

@ctg_names = sort {
	my ( $aComp, $aSubComp, $aSeq ) = ( $a =~ /^comp(\d+)_c(\d+)_seq(\d+)$/ );
	my ( $bComp, $bSubComp, $bSeq ) = ( $b =~ /^comp(\d+)_c(\d+)_seq(\d+)$/ );
	$aComp <=> $bComp || $aSubComp <=> $bSubComp || $aSeq <=> $bSeq;
} @ctg_names;

foreach my $ctg_name ( @ctg_names ) {
	# skip reverse complementary contigs
	#if ( $contig{$ctg_name}->{rc_contig} == 1 ) {
	#	print "$ctg_name\trc_contig:1\n";
	#	next;
	#}

	# for the rest contigs, print all data in a single line per contig
	my ( $comp_name ) = ( $ctg_name =~ /^(comp\d+)_c/ );
	
	print "ctg_name:", $ctg_name, "\t";
	
	print "rc_contig:", $contig{$ctg_name}->{rc_contig}, "\t";
	
	print "length:", $contig{$ctg_name}->{length},  "\t";
	print "length_max:", $comp_max{$comp_name}->{length_max}, "\t";
	
	print "log2ratio:", $contig{$ctg_name}->{log2ratio},  "\t";
	print "log2ratio_max:", $comp_max{$comp_name}->{log2ratio_max}, "\t";
	
	print "seq_cov:", $contig{$ctg_name}->{seq_cov}, "\t";
	
	if ( $contig{$ctg_name}->{seq_cov_dips} ) {
		print "seq_cov_dips:", $contig{$ctg_name}->{seq_cov_dips}, "\t";
		print "seq_cov_dips_max:", $comp_max{$comp_name}->{seq_cov_dips_max}, "\t";
	}
	else {
		# if 'seq_cov_dips' is undef, it means that no reads mapped to it
		# so it should take the maximum penalty
		print "seq_cov_dips:", "1", "\t";
		print "seq_cov_dips_max:", "1", "\t";
	}
	
	print "ctg_proportion_cov_uniref:", $contig{$ctg_name}->{ctg_proportion_cov_uniref}, "\t";
	print "db_proportion_cov_uniref:", $contig{$ctg_name}->{db_proportion_cov_uniref}, "\t";
	print "ctg_proportion_cov_nr:", $contig{$ctg_name}->{ctg_proportion_cov_nr}, "\t";
	print "db_proportion_cov_nr:", $contig{$ctg_name}->{db_proportion_cov_nr}, "\t";
	
	print "strand_uniref:", $contig{$ctg_name}->{strand_uniref}, "\t";
	print "strand_nr:", $contig{$ctg_name}->{strand_nr}, "\t";
	
	print "ctg_length_cov_uniref:", $contig{$ctg_name}->{ctg_length_cov_uniref},  "\t";
	print "ctg_length_cov_uniref_max:", $comp_max{$comp_name}->{ctg_length_cov_uniref_max}, "\t";

	print "db_length_cov_uniref:", $contig{$ctg_name}->{db_length_cov_uniref},  "\t";
	print "db_length_cov_uniref_max:", $comp_max{$comp_name}->{db_length_cov_uniref_max}, "\t";
	
	print "ctg_length_cov_nr:", $contig{$ctg_name}->{ctg_length_cov_nr},  "\t";
	print "ctg_length_cov_nr_max:", $comp_max{$comp_name}->{ctg_length_cov_nr_max}, "\t";

	print "db_length_cov_nr:", $contig{$ctg_name}->{db_length_cov_nr},  "\t";
	print "db_length_cov_nr_max:", $comp_max{$comp_name}->{db_length_cov_nr_max}, "\t";
	
	# Print data related to the called ORFs
	if ( defined ( @{ $contig{$ctg_name}->{orfs} } ) ) {
		my @contig_orfs = @{ $contig{$ctg_name}->{orfs} };
		
		print "orf_number:", scalar(@contig_orfs), "\t";
		
		my @number_of_different_matches = ();
		my $number_of_plus_ORFs = 0;
		my $number_of_minus_ORFs = 0;
		
		my $number_of_plus_ORFs_with_match = 0;
		my $number_of_minus_ORFs_with_match = 0;
		
		
		for ( my $i = 0; $i < @contig_orfs; $i++ ) {
			my @f = split ( /\t/, $contig_orfs[$i] );
			# 0: orf coords
			# 1: orf strand
			# 2: algorithm
			# 3: accession number of the hit
			if ( $f[1] eq '+' ) {
				$number_of_plus_ORFs++;
				if ( $f[3] ne "No_hits" ) {
					$number_of_plus_ORFs_with_match++;
				}
			}
			else {
				$number_of_minus_ORFs++;
				if ( $f[3] ne "No_hits" ) {
					$number_of_minus_ORFs_with_match++;
				}
			}
			if ( $f[3] ne "No_hits" ) {
				push ( @number_of_different_matches, $f[3] );
			}
		}
		
		@number_of_different_matches = uniq(@number_of_different_matches);
		
		print "number_of_plus_ORFs:", $number_of_plus_ORFs, "\t";
		print "number_of_minus_ORFs:", $number_of_minus_ORFs, "\t";
		print "number_of_plus_ORFs_with_match:", $number_of_plus_ORFs_with_match, "\t";
		print "number_of_minus_ORFs_with_match:", $number_of_minus_ORFs_with_match, "\t";
		print "number_of_different_matches:", scalar(@number_of_different_matches);
	}
	else {
		print "orf_number:", "0", "\t";
		print "number_of_plus_ORFs:", "0", "\t";
		print "number_of_minus_ORFs:", "0", "\t";
		print "number_of_plus_ORFs_with_match:", "0", "\t";
		print "number_of_minus_ORFs_with_match:", "0", "\t";
		print "number_of_different_matches:", "0";
	}

	
	print "\n";
}
########################################################################



### Subroutines ########################################################

#load Trinity contigs
sub LoadContigs {
	my $input_file = $_[0];
	open ( TRINITY, $input_file ) or die;
	while ( my $line = <TRINITY> ) {
		if ( $line =~ /^>(comp\d+_c\d+_seq\d+) len=(\d+)/ ) {
			my $ctg_name = $1;
			my $length = $2;
			$contig{$ctg_name} = {
				'length' => $length,
				'rc_contig' => 0,
				'orfs' => [],
				'log2ratio' => 0,
				'seq_cov' => 0,
				'ctg_proportion_cov_uniref' => 0,
				'db_proportion_cov_uniref' => 0,
				'ctg_length_cov_uniref' => 0,
				'db_length_cov_uniref' => 0,
				'strand_uniref' => 0,
				'ctg_proportion_cov_nr' => 0,
				'db_proportion_cov_nr' => 0,
				'ctg_length_cov_nr' => 0,
				'db_length_cov_nr' => 0,
				'strand_nr' => 0,
				'seq_cov_dips' => undef,
			};
			
			# also initialize the hash that will hold certain maximum values
			# per component
			my ( $comp_name ) = ( $ctg_name =~ /^(comp\d+)_/ );
			
			unless ( exists ($comp_max{$comp_name}) ) {
				$comp_max{$comp_name} = {
					'db_length_cov_uniref_max' => 0,
					'ctg_length_cov_uniref_max' => 0,
					'db_length_cov_nr_max' => 0,
					'ctg_length_cov_nr_max' => 0,
					'log2ratio_max' => 0,
					'length_max' => 0,
					'seq_cov_dips_max' => undef,
				};
			}
			
			# check if the length of this contig is the biggest for this component
			if ( $comp_max{$comp_name}->{length_max} ) {
				if ( $contig{$ctg_name}->{length} > $comp_max{$comp_name}->{length_max} ) {
					$comp_max{$comp_name}->{length_max} = $contig{$ctg_name}->{length};
				}
			}
			else {
				$comp_max{$comp_name}->{length_max} = $contig{$ctg_name}->{length};
			}
		}
	}
	close (TRINITY);
}

# load ORFs
sub LoadORFs {
	my $input_file = $_[0];
	open ( ORFS, $input_file ) or die;
	while ( my $line = <ORFS> ) {
		if ( $line =~ /^comp/ ) {
			# skip empty lines
			chomp ($line);
			my @f = split ( /\t/, $line );
			push ( @{ $contig{$f[0]}->{orfs} }, "$f[1]\t$f[2]\t$f[3]\tNo_hits" );
		}
	}
	close (ORFS);
}

# load log2(ratio) of the contigs as well as their sequencing coverage
sub LoadLogRatiosSeqCov {
	my $input_file = $_[0];
	
	open ( RATIOS, $input_file ) or die;
	while ( my $line = <RATIOS> ) {
		chomp ($line);
		
		# add the log2(ratio)
		my @f = split ( /\t/, $line );
		my ( $name, $length ) = ( $f[0] =~/(comp\d+_c\d+_seq\d+)_len_(\d+)$/ );
		
		if ( $f[-1] eq 'N/A' ) {
			$f[-1] = 0;
		}
		
		$contig{$name}->{log2ratio} = $f[-1];
		
		# check if this is the biggest abs(log2ratio) for this component
		my ( $comp_name ) = ( $name=~/^(comp\d+)_c/ );
		if ( $comp_max{$comp_name}->{log2ratio_max} ) {
			if ( abs($f[-1]) > $comp_max{$comp_name}->{log2ratio_max} ) {
				$comp_max{$comp_name}->{log2ratio_max} = abs($f[-1]);
			}
		}
		else {
			$comp_max{$comp_name}->{log2ratio_max} = abs($f[-1]);
		}
		
		# add the sequencing coverage
		$contig{$name}->{seq_cov} = $f[1] * 100 / $length;
	}
	close (RATIOS);
}

# load results from BLASTX against nr and Uniref100
sub LoadBlastxResults {
	my $input_file = $_[0];
	
	my $blast_db = $_[1];
	
	open ( BLASTX, $input_file ) or die;
	
	while ( my $line = <BLASTX> ) {
		chomp ($line);
		my @f = split ( /\t/, $line );
		$contig{$f[0]}->{'ctg_proportion_cov_' . $blast_db} = $f[8]/100;
		$contig{$f[0]}->{'db_proportion_cov_' . $blast_db} = $f[9]/100;
		$contig{$f[0]}->{'ctg_length_cov_' . $blast_db} = $f[6];
		$contig{$f[0]}->{'db_length_cov_' . $blast_db} = $f[7];
		
		if ( $f[-1] eq "Plus" ) {
			$contig{$f[0]}->{'strand_' . $blast_db} = 1;
		}
		else {
			$contig{$f[0]}->{'strand_' . $blast_db} = -1;
		}
		
		# check if "ctg_length_cov" are the biggest for this component
		my ( $comp_name ) = ( $f[0]=~/^(comp\d+)_c/ );
		
		if ( $comp_max{$comp_name}->{'ctg_length_cov_' . $blast_db . '_max'} ) {
			if ( $f[6] > $comp_max{$comp_name}->{'ctg_length_cov_' . $blast_db . '_max'} ) {
				$comp_max{$comp_name}->{'ctg_length_cov_' . $blast_db . '_max'} = $f[6];
			}
		}
		else {
			$comp_max{$comp_name}->{'ctg_length_cov_' . $blast_db . '_max'} = $f[6];
		}
		
		# same for "db_length_cov"
		if ( $comp_max{$comp_name}->{'db_length_cov_' . $blast_db . '_max'} ) {
			if ( $f[7] > $comp_max{$comp_name}->{'db_length_cov_' . $blast_db . '_max'} ) {
				$comp_max{$comp_name}->{'db_length_cov_' . $blast_db . '_max'} = $f[7];
			}
		}
		else {
			$comp_max{$comp_name}->{'db_length_cov_' . $blast_db . '_max'} = $f[7];
		}
	}
	
	close (BLASTX);
}

sub LoadOrfBlastpResults {
	# not all ORFs present in the BLAST results should be present in the contigs
	my $input_file = $_[0];
	my $algorithm = $_[1];
	
	open ( BLASTP, $input_file ) or die;
	
	while ( my $line = <BLASTP> ) {
		chomp ($line);
		my @f = split ( /\t/, $line );
		
		my ( $ctg_name, $orf_coords, $orf_strand ) = ( $f[0] =~ /^(comp\d+_c\d+_seq\d+):(\d+-\d+)\(([-\+])\)$/ );
		
		# look if this ORF is found in the particular contig
		if ( defined ( @{ $contig{$ctg_name}->{orfs} } ) ) {
			my @contig_orfs = @{ $contig{$ctg_name}->{orfs} };
			
			for ( my $i = 0; $i < @contig_orfs; $i++ ) {
				my @ff = split ( /\t/, $contig_orfs[$i] );
				if ( $ff[0] eq $orf_coords && $ff[1] eq $orf_strand && $ff[2] eq $algorithm ) {
					# if you find exactly this ORF in the list of called ORFs
					# for this particular contig, add the accession number of the best hit
					$ff[3] = $f[2];
					$contig_orfs[$i] = join ( "\t", @ff );
				}
			}
			
			@{ $contig{$ctg_name}->{orfs} } = @contig_orfs;
		}
	}
	close (BLASTP);
}

sub TagRCContigs {
	my $input_file = $_[0];
	open ( RC_CONTIGS, $input_file ) or die;
	while ( my $ctg_name = <RC_CONTIGS> ) {
		chomp ($ctg_name);
		$contig{$ctg_name}->{rc_contig} = 1;
	}
	close (RC_CONTIGS);
}

sub LoadSeqCovDips {
	my $input_file = $_[0];
	
	open ( SEQ_COV_DIPS, $input_file ) or die;
	while  ( my $line = <SEQ_COV_DIPS> ) {
		chomp ($line);
		my @f = split ( /\t/, $line );
		my @node_cov = split ( /, /, $f[1] );
		
		my @coverage_dips = ();
		
		for ( my $i = 0; $i < (scalar(@node_cov) - 1); $i++ ) {
			my ( $node1, $node2 ) = ( $node_cov[$i], $node_cov[$i+1] );
			
			my @node1 = split ( /:/, $node1 );
			my @node2 = split ( /:/, $node2 );
			# [0] is node name
			# [1] is node coords
			# [2] is coverage
			
			if ( $node1[2] == 0 && $node2[2] == 0 ) {
				next; # no change in coverage
			}
			
			if ( $node1[2] < $node2[2] ) {
				if ( $node1[2] / $node2[2] < 0.2 ) {
					push ( @coverage_dips, "$node1 => $node2" );
				}
			}
			else {
				if ( $node2[2] / $node1[2] < 0.2 ) {
					push ( @coverage_dips, "$node1 => $node2" );
				}
			}
		}
		
		$contig{$f[0]}->{'seq_cov_dips'} = scalar (@coverage_dips);
		
		# check if this number of dips is the biggest for this component
		my ( $comp_name ) = ( $f[0]=~/^(comp\d+)_c/ );
		
		if ( $comp_max{$comp_name}->{seq_cov_dips_max} ) {
			if ( scalar (@coverage_dips) > $comp_max{$comp_name}->{seq_cov_dips_max} ) {
				$comp_max{$comp_name}->{seq_cov_dips_max} = scalar (@coverage_dips);
			}
		}
		else {
			$comp_max{$comp_name}->{seq_cov_dips_max} = scalar (@coverage_dips);
		}

	}
	
	close (SEQ_COV_DIPS);
}
