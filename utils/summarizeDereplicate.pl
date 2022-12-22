#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(:config auto_abbrev no_ignore_case);
use FileHandle;
$|=1;

sub usage {
	my $message = shift;
	my $usage = qq/
usage: summarizeDereplicate.pl [options]

Required:
-r, --regions       FILE   regions.txt file (output by starfish dereplicate)
-d, --dereplicated  FILE   dereplicated.txt file (output by starfish dereplicate)

Required, with defaults:
-s, --separator     STR    character separating genomeID from featureID (default: '_')

Optional:
-h, --help              print more details and exit
/;
	if (not defined $message) {
		$message = qq/
Summarizes statistics of dereplicated regions. Useful if .regions.txt and .dereplicate.txt
files had to be manually edited. Prints to STDOUT. \n\n/;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($usage, $message);
}

main : {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'dereplicated|d=s',
		'regions|r=s',
		'separator|s=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	# parse the number of reference elements per group in each region
	# dereplicatedValidatedCount{$regionID}{$groupID} = # of refs with validated sites
	my ($dereplicatedValidatedCount) = Parse_ref_per_region($opts{'dereplicated'});
	
	# parse regions stats
	my ($region2groupCount, $region2emptyCount, $region2fragmentedCount, $region2elementCount) = Parse_region_data($opts{'regions'});
	
	# $region2groupCount{$regionID}{$group}{elementID} = 1
	# region2emptyCount{$regionID}{contigID} = 1
	# region2fragmentedCount{$regionID}{contigID} = 1
	# region2elementCount{$regionID} = # of elements in region
	Print_summary_report($region2groupCount, $region2emptyCount, $region2fragmentedCount, $dereplicatedValidatedCount, $region2elementCount, $opts{'separator'});


}

sub Print_summary_report {
	my ($region2elementGroup, $region2empty, $region2fragmented, $dereplicatedValidatedCount, $region2elementCount, $SEP) = @_;
	print "#regionID\tstarshipFamily\tstarshipElements\trefsWithValidatedSites\temptyRegions\tfragmentedRegions\tomesWithElementGroup\tomesWithOtherElements\tomesWithEmpty\tomesWithFrag\n";

	foreach my $regionID (sort keys %{$region2elementGroup}) {
		
		my (%genomeElementCount, %genomeEmptyCount, %genomeFragCount); # keep track of how many genomes have at least 1 of each allele at this region
		my ($elementRegionCount, $emptyRegionCount, $fragmentedRegionCount) = (0,0,0);

		# gather fragmented info
		foreach my $fragmentedContigID (sort keys %{$region2fragmented->{$regionID}}) {
			$fragmentedRegionCount++;
			my ($omeID) = split/$SEP/, $fragmentedContigID;
			$genomeFragCount{$omeID} = 1;
		}
		
		# gather empty info (ignore any genomes with empty that also have a frag; old bug from starfish dereplicate)
		foreach my $emptyContigID (sort keys %{$region2empty->{$regionID}}) {
			my ($omeID) = split/$SEP/, $emptyContigID;
			if (not exists $genomeFragCount{$omeID}) {
				$emptyRegionCount++;
				$genomeEmptyCount{$omeID} = 1;
			}
		}

		# figure out how many elements and genomes with elements total
		foreach my $elementID (sort keys %{$region2elementCount->{$regionID}}) {
			$elementRegionCount++;
			my ($omeID) = split/$SEP/, $elementID;
			$genomeElementCount{$omeID} = 1;
		}

		my $omesWithEmpty = scalar (keys %genomeEmptyCount);
		my $omesWithFrag = scalar (keys %genomeFragCount);

		# gather element info, and print out group-specific info
		foreach my $elementGroupID (sort keys %{$region2elementGroup->{$regionID}}) {
			my %genomeElementGroupCount;
			my $elementGroupCount = 0;
			foreach my $elementID (sort keys %{$region2elementGroup->{$regionID}->{$elementGroupID}}) {
				$elementGroupCount++;
				my ($omeID) = split/$SEP/, $elementID;
				$genomeElementGroupCount{$omeID} = 1;
			}
			my $validatedCount = 0;
			$validatedCount = $dereplicatedValidatedCount->{$regionID}->{$elementGroupID} if (exists $dereplicatedValidatedCount->{$regionID}->{$elementGroupID});
			my $omesWithElement = scalar (keys %genomeElementGroupCount);
			my $omesWithOtherElements = scalar(keys %genomeElementCount) - $omesWithElement;
			print "$regionID\t$elementGroupID\t$elementGroupCount\t$validatedCount\t$emptyRegionCount\t$fragmentedRegionCount\t$omesWithElement\t$omesWithOtherElements\t$omesWithEmpty\t$omesWithFrag\n";		
		}
	}	
}


sub Parse_ref_per_region {
	my ($dereplicatedFile) = @_;
	
	# parse reference element
	# if a region has multiple reference elements, pick the one associated with the most elements
	my (%validatedRefs);
	open(my $derepIN, '<', $dereplicatedFile) or usage("\nError: can't open $dereplicatedFile for reading\n");
	while (my $line = <$derepIN>) {
		chomp $line;
		next if ($line =~ m/^#/);
		my ($regionID, $elementGroupID, $refElementID, $refSiteID, $otherElementIDs) = split/\t/, $line;
		if ($refSiteID ne '.') {
			$validatedRefs{$regionID}{$elementGroupID}++;
		}
	}
	return(\%validatedRefs);
}

sub Parse_region_data {
	my ($regionFile) = @_;

	# parse region data
	my (%region2groupCount, %region2emptyCount, %region2fragmentedCount, %region2elementCount);
	open(my $regionIN, '<', $regionFile) or usage("\nError: can't open $regionFile for reading\n");
	while (my $line = <$regionIN>) {
		chomp $line;
		next if ($line =~ m/^#/);
		my ($regionID, $memberGroupID, $memberID, $memberType, $contigID, $begin, $end, $regionBegin, $regionEnd, $flankingOGs) = split/\t/, $line;
		
		if ($memberType eq 'frag') { 
			$region2fragmentedCount{$regionID}{$contigID} = 1;
		} elsif ($memberType eq 'empty') { 
			$region2emptyCount{$regionID}{$contigID} = 1;
		} else {
			$region2groupCount{$regionID}{$memberGroupID}{$memberID} = 1;
			$region2elementCount{$regionID}{$memberID} = 1;
		}
	}
	return(\%region2groupCount, \%region2emptyCount, \%region2fragmentedCount, \%region2elementCount);
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a file to --dereplicated\n") if (not defined $opts->{'dereplicated'});
	usage("\nError: the file provided to --dereplicated does not exist\n") if (! -f $opts->{'dereplicated'});
	usage("\nError: please provide a file to --regions\n") if (not defined $opts->{'regions'});
	usage("\nError: the file provided to --regions does not exist\n") if (! -f $opts->{'regions'});
	if (not defined $opts->{'separator'}) {
		$opts->{'separator'} = '_';
	} elsif ($opts->{'separator'} eq ':') {
		usage("\nError: the separator character cannot be ':'\n");
	} elsif ($opts->{'separator'} eq ';') {
		usage("\nError: the separator character cannot be ';'\n");
	} elsif ($opts->{'separator'} eq '|') {
		usage("\nError: the separator character cannot be '|'\n");
	}
	if (defined $opts->{'separator'}) {
		$opts->{'separator'} = quotemeta($opts->{'separator'});	# to allow splitting on special characters, like '.'
	}
}



