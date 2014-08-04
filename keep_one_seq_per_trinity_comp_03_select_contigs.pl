#!/usr/local/bin/perl -w

use warnings;
use strict;
#use List::MoreUtils qw(uniq);

if( scalar ( @ARGV ) != 2 ){
	my @temp=split(/\//, $0);
	die
"

USAGE: $temp[-1] <List of contigs with their scores> <score cutoff>

This script reads the output file of the previous script ('*calculate_scores.pl')
and selects at most one contig per component.

";
}

my $SCORE_CUTOFF = $ARGV[1];

my %component = ();

open ( CONTIG_SCORES, $ARGV[0] ) or die;
while ( my $line = <CONTIG_SCORES> ) {
	chomp ($line);
	my @f = split ( /\t/, $line );
	
	# discard low-scoring contigs
	if ( $f[2] <= $SCORE_CUTOFF ) { next; }
	
	my ( $comp_name ) = ( $f[0] =~ /^(comp\d+)_/ );
	if ( !exists ( $component{$comp_name} ) ) {
		$component{$comp_name} = [];
	}
	
	push ( @{ $component{$comp_name} }, $line );
	
}
close (CONTIG_SCORES);

#print scalar ( keys (%component) ), "\n";

my @sorted_component_names = sort {
	#print STDERR "$a <=> $b\n";
	#my $asdf = <STDIN>;
	my ( $aVal ) = ( $a =~ /^comp(\d+)$/ );
	my ( $bVal ) = ( $b =~ /^comp(\d+)$/ );
	$aVal <=> $bVal;
} ( keys ( %component ) );

foreach my $comp_name ( @sorted_component_names ) {
	my @comp_contigs = @{ $component{$comp_name} };
	#print STDERR $comp_contigs[0], "\n";
	#my $asdf = <STDIN>;
	
	my @sorted_comp_contigs = sort {
		my @fa = split ( /\t/, $a );
		my @fb = split ( /\t/, $b );
		
		$fb[2] <=> $fa[2];
	} @comp_contigs;
	
	print $sorted_comp_contigs[0], "\n";
}
