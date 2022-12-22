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
usage: neighborhoods2mat.pl [options]

Required:
-b, --bed     FILE   BED file with neighborhoods (e.g., output by glofish sketch or starfish summarize)
-t, --tag     FILE   1 column tsv with: tag of interest (corresponding to tag column in BED)

Optional:
-g, --gene    FILE   1 column tsv with: focal tag(s) for printing gene x tag matrix
-h, --help           print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script reformats a BED file into matrix format and prints a region x tag matrix by
default. If --gene file provided, will print a gene x tag matrix instead of a region x
tag matrix. Prints to STDOUT. \n\n/;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($usage, $message);
}

main: {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'bed|b=s',
		'tag|t=s',
		'gene|g=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	# load up annotations of interest
	my ($focalTags) = dim_0_hash($opts{'tag'}, "\t", "0");
	my ($geneTags) = dim_0_hash($opts{'gene'}, "\t", "0") if (defined $opts{'gene'});

	# format data into matrix format
	# {hoodID/focalGeneID}{geneBegin}{geneID} = tag
	my ($neighborhoods) = Parse_bed($opts{'bed'}, $geneTags);
	
	# print a neighborhood/gene x tag genotyping matrix
	Print_matrix($neighborhoods, $focalTags);
}

sub Print_matrix {
	my ($queryNeighborhoods, $focalTags) = @_;
	
	# print header
	print "#Names";
	foreach my $tag (sort keys %{$focalTags}) {
		print "\t$tag";
	}
	print "\n";
	
	# count up all genes with annotations of interest per neighborhood
	my %hoodAnnotations;
	foreach my $hoodID (keys %{$queryNeighborhoods}) {
		foreach my $geneBegin (keys %{$queryNeighborhoods->{$hoodID}}) {
			foreach my $geneID (keys %{$queryNeighborhoods->{$hoodID}->{$geneBegin}}) {
				my $queryTag = $queryNeighborhoods->{$hoodID}->{$geneBegin}->{$geneID};
				if (exists $focalTags->{$queryTag}) {
					$hoodAnnotations{$hoodID}{$queryTag}++;
				}
			}
		}
	}

	# print out occupancy matrix
	foreach my $hoodID (sort keys %hoodAnnotations) {
		print "$hoodID";
		foreach my $tag (sort keys %{$focalTags}) {
			if (exists $hoodAnnotations{$hoodID}{$tag}) {
				print "\t$hoodAnnotations{$hoodID}{$tag}";
			} else {
				print "\t0";
			}
		}
		print "\n";
	}
}

sub Parse_bed {
	my ($neighborhoodFile, $geneTags) = @_;
	#   structured:  $neighborhoods{$neighborhoodID}{$begin}{$queryID} = $tag
	my (%neighborhood2focalGene, %neighborhoods);
	
	open (my $IN, '<', $neighborhoodFile) or usage("\nError: could not open $neighborhoodFile for reading\n");
	while (my $line = <$IN>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($contigID, $begin,$end, $featID, $featTag, $featStrand, $neighborhoodID) = split/\t/, $line;
		$neighborhoods{$neighborhoodID}{$begin}{$featID} = $featTag;
		if (defined $geneTags && exists $geneTags->{$featTag}) {
			$neighborhood2focalGene{$neighborhoodID}{$featID} = 1;
		}
	}
	if (defined $geneTags) {
		my %geneNeighborhood;
		foreach my $hoodID (keys %neighborhood2focalGene) {
			foreach my $focalGeneID (keys %{$neighborhood2focalGene{$hoodID}}) {
				foreach my $geneBegin (keys %{$neighborhoods{$hoodID}}) {
					foreach my $geneID (keys %{$neighborhoods{$hoodID}{$geneBegin}}) {
						my $queryTag = $neighborhoods{$hoodID}{$geneBegin}{$geneID};
						$geneNeighborhood{$focalGeneID}{$geneBegin}{$geneID} = $queryTag;
					}
				}
			}
		}
		return(\%geneNeighborhood);
	
	} else {
		return(\%neighborhoods);	
	}
}

sub dim_0_hash{
my $usage = qq/
usage: my (\%hash) =  dim_0_hash(<abs_path_to_file>, <separator>, <column_numbers>)
sub 0_dim_hash inputs the path to a tab delimited file, and returns a hash where key = column number and value = 1\n\n/;
	my ($file, $sep, $coi) = @_;
	die($usage, "Error: incorrect number of args..\n") if (scalar @_ != 3);
	die($usage, "Error: only one column number should be specified\n") if ($coi =~ m/[^\d]/);
	my %hash;
	open(my $in, '<', $file) or die($usage, "Error: cannot open $file..\n");
	while (my $line = <$in>) {
		if ($line =~ m/^#/) {
			next;
		} else {
			chomp $line;
			my @cols = split/$sep/, $line;
			my ($first) = ($cols[$coi]);
			$hash{$first} = 1;
		}
	}
	return(\%hash);
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a file to --bed\n") if (not defined $opts->{'bed'});
	usage("\nError: the file provided to --bed does not exist\n") if (! -f $opts->{'bed'});
	usage("\nError: please provide a file to --tag\n") if (not defined $opts->{'tag'});
	usage("\nError: the file provided to --tag does not exist\n") if (! -f $opts->{'tag'});
	if (defined $opts->{'gene'}) {
		usage("\nError: the file provided to --gene does not exist\n") if (! -f $opts->{'gene'});
	}
}