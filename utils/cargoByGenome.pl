#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Std;
use FileHandle;
$|=1;

sub usage {
	my $message = shift;
	my $usage = qq/
usage: $0 [OPTIONS]

Required:
-o    FILE    orthofinder-formatted orthogroups.txt file of all ortholog groups in analysis
-g    FILE    2 column tsv with: genomeID, category of interest (e.g., plant host range)
-e    FILE    starships.bed file (output by starfish summarize)

Required, with defaults:
-s    STR     character separating genomeID from featureID (default: '_')\n/;

	if (not defined $message) {
		$message = qq/
This script prints to STDOUT out a long-form summary of the number of genomes per genome
category with Starship-associated ortholog groups\n\n/; 
	} else {
		$message = "\n$message\nuse -h for more details\n" ;
	}	
	die($usage, $message);
}

main: {
	my %opts;
	getopt('o:g:e:s:h', \%opts);
	Opts_check(\%opts);
	
	# read in data
	my ($ome2cat) = dim_1_hash($opts{'g'}, "\t", "0:1");
	my ($og2gene) = Parse_group_file_by_og($opts{'o'});
	my ($gene2og) = Parse_group_file_by_member($opts{'o'});
	
	# grab all OGs in starships
	my $starshipGenes = dim_0_hash($opts{'e'}, "\t", "3");
	my %starshipOGs;
	foreach my $geneID (keys %{$starshipGenes}) {
		$starshipOGs{$gene2og->{$geneID}} = 1 if (exists $gene2og->{$geneID});
	}
	
	# count number of starship-associated OGs in starships and not in starships for each genome category
	my $datestring = localtime();					
	my %ogCounts;
	foreach my $og (keys %starshipOGs) {
		foreach my $geneID (@{$og2gene->{$og}}) {
			my ($omeID) = split/$opts{s}/, $geneID;
			
			if (exists $ome2cat->{$omeID}) {
				if (exists $starshipGenes->{$geneID}) {
					$ogCounts{$og}{$ome2cat->{$omeID}}{$omeID}{'starship'} = 1;
				} else {
					$ogCounts{$og}{$ome2cat->{$omeID}}{$omeID}{'background'} = 1;
				}
			} else {
				usage("\n\n[$datestring] error: could not find a category for genome $omeID, exiting\n");
			}
		}
	}
	
	# print, account for genomes that may have both starship and genome associated members of the same OG
	print "ogID\tcategory\tcompartment\tgenomeCount\n";
	foreach my $og (sort keys %ogCounts) {
		foreach my $category (sort keys %{$ogCounts{$og}}) {
			my ($starshipCount, $backgroundCount, $bothCount) = (0,0,0);
			foreach my $omeID (keys %{$ogCounts{$og}{$category}}) {
				if (exists $ogCounts{$og}{$category}{$omeID}{'starship'} && exists $ogCounts{$og}{$category}{$omeID}{'background'}) {
					$bothCount++;
				} elsif (exists $ogCounts{$og}{$category}{$omeID}{'starship'}) {
					$starshipCount++;
				} elsif (exists $ogCounts{$og}{$category}{$omeID}{'background'}) {
					$backgroundCount++;
				}
			}
			print "$og\t$category\tstarship\t$starshipCount\n";
			print "$og\t$category\tbackground\t$backgroundCount\n";
			print "$og\t$category\tboth\t$bothCount\n";
		}
	}	
}

sub Parse_group_file_by_member {
	my ($clusteringOutfile) = @_;
	my $datestring = localtime();					
	my %member2group;
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cannot read $clusteringOutfile, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		my (@members) = split/\s+/, $line;
		my $group = shift @members;
		$group =~ s/:$//;
		
		foreach my $member (@members) {
			$member2group{$member} = $group;
		}
	}
	return(\%member2group);	
}

sub Parse_group_file_by_og {
	my ($clusteringOutfile) = @_;
	my $datestring = localtime();					
	my (%group2member);
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cant open $clusteringOutfile for reading, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		my (@members) = split/\s+/, $line;
		my $group = shift @members;
		$group =~ s/:$//;
		push @{$group2member{$group}}, @members;
	}
	return(\%group2member);	
}


sub dim_1_hash{
	my ($file, $sep, $columnstring) = @_;
	my @coi = split/:/, $columnstring; #columns of interest
	my %hash;
	my $datestring = localtime();					
	open(my $in, '<', $file) or usage("\n\n[$datestring] error: cant open $file for reading, exiting\n");
	while (my $line = <$in>) {
		if ($line =~ m/^#/) {
			next;
		} else {
			chomp $line;
			my @cols = split/$sep/, $line;
			my ($first, $second) = ($cols[$coi[0]], $cols[$coi[1]]);
			usage("Error: cannot parse line from $file:\n$line\n") if (not defined $first);
			usage("Error: cannot parse line from $file:\n$line\n") if (not defined $second);
			$hash{$first} = $second;
		}
	}
	return(\%hash);
}

sub dim_0_hash{
	my ($file, $sep, $coi) = @_;
	my %hash;
	open(my $in, '<', $file) or usage("Error: can't open $file for reading\n");
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
	usage() if (defined $opts->{'h'});
	usage("\nError: please provide a file to -o\n") if (not defined $opts->{'o'});
	usage("\nError: the file provided to -o does not exist\n") if (! -f $opts->{'o'});
	usage("\nError: please provide a file to -g\n") if (not defined $opts->{'g'});
	usage("\nError: the file provided to -g does not exist\n") if (! -f $opts->{'g'});
	usage("\nError: please provide a file to -e\n") if (not defined $opts->{'e'});
	usage("\nError: the file provided to -e does not exist\n") if (! -f $opts->{'e'});
	if (not defined $opts->{'s'}) {
		$opts->{'s'} = '_';
	} elsif ($opts->{'s'} eq ':') {
		usage("\nError: the s character cannot be ':'\n");
	} elsif ($opts->{'s'} eq ';') {
		usage("\nError: the s character cannot be ';'\n");
	} elsif ($opts->{'s'} eq '|') {
		usage("\nError: the s character cannot be '|'\n");
	}
	if (defined $opts->{'s'}) {
		$opts->{'s'} = quotemeta($opts->{'s'});	# to allow splitting on special characters, like '.'
	}
}
