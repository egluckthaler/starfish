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
usage: regions2mat.pl [options]

Required:
-r, --regions     FILE   .regions.txt file with region metadata (output by starfish dereplicate)
-m, --mode        STR    type of matrix to print, either 'starship' or 'genome'

Required, if --mode genome
-s, --separator   STR    character separating genomeID from featureID

Optional:
-h, --help               print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script reformats the .regions.txt file output by starfish dereplicate into matrix 
format. Mostly useful for manipulating subsetted .regions.txt files containing regions
of interest. Print to STDOUT\n\n/;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($usage, $message);
}

main: {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'regions|r=s',
		'mode|m=s',
		'separator|s=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	# format data into matrix
	my ($starship2region, $genome2region, $regionIDs) = Format_region_matrix($opts{'regions'}, $opts{'separator'});

	if ($opts{'mode'} eq 'starship') {
		# print a starship x region genotyping matrix
		Print_matrix($starship2region, $regionIDs);
	} elsif ($opts{'mode'} eq 'genome') {
		# print a starship x region genotyping matrix
		Print_matrix($genome2region, $regionIDs);
	}	
}

sub Print_matrix {
	my ($row2region, $regionIDs, $outfile) = @_;
	
	# print header
	print "#Names";
	foreach my $regionID (sort keys %{$regionIDs}) {
		print "\t$regionID";
	}
	print "\n";
	
	# print row data
	foreach my $rowID (sort keys %{$row2region}) {
		print "$rowID";
		foreach my $regionID (sort keys %{$regionIDs}) { # make sure every region gets printed
			if (exists $row2region->{$rowID}->{$regionID}) {
				print "\t$row2region->{$rowID}->{$regionID}";
			} else {
				print "\t0";
			}
		}
		print "\n";
	}
}

sub Format_region_matrix {
	my ($regionFile, $SEP) = @_;

	my (%starship2region, %ome2region, %regionIDs);
	
	open (my $IN, '<', $regionFile) or usage("\nError: could not open $regionFile for reading\n");
	
	while (my $line = <$IN>) {
		next if ($line =~ m/^#/);
		chomp $line;
		
		my ($regionID,$memberGroupID,$memberID,$contigID,$begin,$end,$regionBegin,$regionEnd,$flankingOGs) = split/\t/, $line;
		
		$regionIDs{$regionID} = 1;
		
		if ($memberGroupID ne '.') { # only starships are assigned to families
	
			# element info, where element symbolized by 1
			$starship2region{$memberID}{$regionID} = 1;
			
			if (defined $SEP) {
				my ($omeID)= split/$SEP/, $contigID;
				$ome2region{$omeID}{$regionID} = 2;
			}
		} elsif ($flankingOGs !~ m/,\.,/) { # only empty regions/sites dont have a . marker in OGs; 

			# empty info, where empty symbolized by 0
			if (defined $SEP) {
				my ($omeID)= split/$SEP/, $contigID;
				$ome2region{$omeID}{$regionID} = -1;
			}
		} else { # everything else is a fragment

			if (defined $SEP) {
				my ($omeID)= split/$SEP/, $contigID;
				$ome2region{$omeID}{$regionID} = 1;
			}
		}
	}
	return(\%starship2region, \%ome2region, \%regionIDs);	
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a .regions.txt file for --regins\n") if (not defined $opts->{'regions'});
	usage("\nError: the .regions.txt file provided to --regions does not exist\n") if (! -f $opts->{'regions'});
	usage("\nError: please specify an option for --mode\n")  if (not defined $opts->{'mode'});
	usage("\nError: you must specify either 'starship' or 'genome' as an option for --mode\n")  if ($opts->{'mode'} !~ m/^genome$|^starship$/);
	if ($opts->{'mode'} eq 'genome') {
		if (defined $opts->{'separator'}) {
			if ($opts->{'separator'} eq ':') {
				usage("\nError: the separator character cannot be ':'\n");
			} elsif ($opts->{'separator'} eq ';') {
				usage("\nError: the separator character cannot be ';'\n");
			} elsif ($opts->{'separator'} eq '|') {
				usage("\nError: the separator character cannot be '|'\n");
			}
			$opts->{'separator'} = quotemeta($opts->{'separator'});	# to allow splitting on special characters, like '.'
		} else {
			usage("\nError: you must provide a --separator if using --mode genome\n");
		}
	}
}