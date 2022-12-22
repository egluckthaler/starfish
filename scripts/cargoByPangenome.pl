#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Std;
use FileHandle;
$|=1;


# change global vars here
my $MINREPEAT = 200;

sub usage {
	my $message = shift;
	my $usage = qq/
usage: $0 [OPTIONS]

Required:
-o    FILE    orthofinder-formatted orthogroups.txt file of all ortholog groups in analysis
-b    FILE    starships.bed file (output by starfish summarize)
-g    INT     total number of genomes in analysis

Required, with defaults:
-s    STR     character separating genomeID from featureID (default: '_')

Optional:
-e    FILE    1 column tsv with genome codes to ignore
-r    FILE    2 column tsv with: genome code, path to GFF of repetitive sequence coordinates\n/;

	if (not defined $message) {
		$message = qq/
This script prints to STDOUT out a long-form summary of the number of genomes with and 
without each OG, and the number of genomes that OG is found in Starships, genomic 
background, or both. Any ortholog group where at least 1 member overlaps with a repetitive
sequence will be removed\n\n/; 
	} else {
		$message = "\n$message\nuse -h for more details\n" ;
	}	
	die($usage, $message);
}

main: {
	my %opts;
	getopt('o:b:s:g:r:e:h', \%opts);
	Opts_check(\%opts);
	
	# read in data
	my ($omesToIgnore) = dim_0_hash($opts{'e'}, "\t", "0");
	my ($og2gene) = Parse_group_file_by_og($opts{'o'}, $omesToIgnore);
	my ($gene2og) = Parse_group_file_by_member($opts{'o'}, $omesToIgnore);

	my ($restrictedCoords);
	if (defined $opts{'r'}) {
		my ($ome2repeatGFF) = dim_1_hash($opts{'r'}, "\t", "0:1");
		($restrictedCoords) = Parse_repeat_GFF($ome2repeatGFF);
	}


	# Parse genes in starships
	# structured: {geneID}{contigID} = (begin, end)
	my ($starshipGenes) = Parse_feature_coords($opts{'b'});
	
	# identify OGs to ignore
	my %restrictedOGs;
	foreach my $geneID (keys %{$starshipGenes}) {
		foreach my $contigID (keys %{$starshipGenes->{$geneID}}) {
			my ($featBegin, $featEnd) = @{$starshipGenes->{$geneID}->{$contigID}};
			my $restrictedFound = 0;
			if (defined $restrictedCoords) {
				if (exists $restrictedCoords->{$contigID}) {
					foreach my $restrictedCoord (keys %{$restrictedCoords->{$contigID}}) {
						if ((($restrictedCoord >= $featBegin) && ($restrictedCoord <= $featEnd))) {
							$restrictedFound = 1;
							last;
						}
					}
				}
	
			}
			# identify OGs to ignore
			if ($restrictedFound == 1) {
				$restrictedOGs{$gene2og->{$geneID}} = 1 if (exists $gene2og->{$geneID});
			}
		}
	}
	
	# filter out any OG without a member that overlaps a repeat
	my %filteredOGs;
	foreach my $og (keys %{$og2gene}) {
		$filteredOGs{$og} = 1 if (not exists $restrictedOGs{$og});
	}
	
	# count number of genomes with OGs in starships, in genomic background, in both for each OG
	my $datestring = localtime();					
	my %ogCounts;
	foreach my $og (keys %filteredOGs) { # only look at filtered OGs
		foreach my $geneID (@{$og2gene->{$og}}) {
			my ($omeID) = split/$opts{s}/, $geneID;
			if (exists $starshipGenes->{$geneID}) {
				$ogCounts{$og}{$omeID}{'starship'}++;
			} else {
				$ogCounts{$og}{$omeID}{'background'}++;
			}
		}
	}
	
	# print, account for genomes that may have both starship and genome associated members of the same OG
	print "ogID\tbackgroundGenomePresence\tstarshipGenomePresence\tbothGenomePresence\tcompartmentType\ttotalPresence\ttotalAbsence\tprevalence\tbackgroundMemberCount\tstarshipMemberCount\ttotalCount\tavgCount\tstarshipMemberFraction\tquadrant\n";
	foreach my $og (sort keys %ogCounts) {
		my ($starshipMemberCount, $backgroundMemberCount, $starshipCount, $backgroundCount, $bothCount) = (0,0,0,0,0);
		foreach my $omeID (keys %{$ogCounts{$og}}) {
			if (exists $ogCounts{$og}{$omeID}{'starship'} && exists $ogCounts{$og}{$omeID}{'background'}) {
				$bothCount++;
				$starshipMemberCount += $ogCounts{$og}{$omeID}{'starship'};
				$backgroundMemberCount += $ogCounts{$og}{$omeID}{'background'};
			} elsif (exists $ogCounts{$og}{$omeID}{'starship'}) {
				$starshipCount++;
				$starshipMemberCount += $ogCounts{$og}{$omeID}{'starship'};
			} elsif (exists $ogCounts{$og}{$omeID}{'background'}) {
				$backgroundCount++;
				$backgroundMemberCount += $ogCounts{$og}{$omeID}{'background'};
			}
		}
		my $totalPresence = $bothCount + $starshipCount + $backgroundCount;
		my $totalAbsence = $opts{'g'} - $totalPresence;
		
		# calculate the minor presence/absence frequency
		my $MAF = sprintf("%.3f", $totalPresence / $opts{'g'});
		$MAF = sprintf("%.3f",  $totalAbsence / $opts{'g'}) if (($totalAbsence / $opts{'g'}) < $MAF);
		
		# print out compartment type:
		my $type = 'background_only';
		if (($starshipCount > 0) || ($bothCount > 0)) {
			$type = 'starship_associated';
		} elsif (($starshipCount > 0) && ($bothCount == 0) && ($backgroundCount == 0)) {
			$type = 'starship_only';
		} 	
		
		# calculate fraction of members found in starships
		my $starshipMemberFraction = sprintf("%.3f", $starshipMemberCount / ($starshipMemberCount + $backgroundMemberCount));

		# calculate OG prevalence pangenome-wide
		my $prevalence = sprintf("%.3f", $totalPresence / ($totalPresence + $totalAbsence));

		# calculate quadrant based on 0.1, 0.9 thresholds grouping
		my ($lowThresh, $highThresh) = (0.1, 0.9);
		
		my $quadrant;
		if (($starshipMemberFraction < $lowThresh) && ($prevalence >= $highThresh)) {
			$quadrant = 1;
		} elsif (($starshipMemberFraction >= $lowThresh) && ($starshipMemberFraction < $highThresh) && ($prevalence >= $highThresh)) {
			$quadrant = 2;
		} elsif (($starshipMemberFraction >= $highThresh) && ($prevalence >= $highThresh)) {
			$quadrant = 3;
		} elsif (($starshipMemberFraction < $lowThresh) && ($prevalence >= $lowThresh) && ($prevalence < $highThresh)) {
			$quadrant = 4;
		} elsif (($starshipMemberFraction >= $lowThresh) && ($starshipMemberFraction < $highThresh) && ($prevalence >= $lowThresh) && ($prevalence < $highThresh)) {
			$quadrant = 5;
		} elsif (($starshipMemberFraction >= $highThresh) && ($prevalence >= $lowThresh) && ($prevalence < $highThresh)) {
			$quadrant = 6;
		} elsif (($starshipMemberFraction < $lowThresh) && ($prevalence < $lowThresh)) {
			$quadrant = 7;
		} elsif (($starshipMemberFraction >= $lowThresh) && ($starshipMemberFraction < $highThresh) && ($prevalence < $lowThresh)) {
			$quadrant = 8;
		} elsif (($starshipMemberFraction >= $highThresh) && ($prevalence < $lowThresh)) {
			$quadrant = 9;
		}
		
		my $totalCount = $backgroundMemberCount + $starshipMemberCount;
		my $avgCount = sprintf("%.3f", $totalCount / $opts{'g'});
		print "$og\t$backgroundCount\t$starshipCount\t$bothCount\t$type\t$totalPresence\t$totalAbsence\t$prevalence\t$backgroundMemberCount\t$starshipMemberCount\t$totalCount\t$avgCount\t$starshipMemberFraction\tq$quadrant\n";
	}	
}

sub Parse_group_file_by_member {
	my ($clusteringOutfile, $excludedOmes) = @_;
	my $datestring = localtime();					
	my %member2group;
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cannot read $clusteringOutfile, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		next if ($line =~ m/^Orthogroup/);
		my (@members) = split/\s+/, $line;
		my $group = shift @members;
		$group =~ s/:$//;
		
		foreach my $member (@members) {
			my $memberFail = 0;
			if (defined $excludedOmes) {
				foreach my $excludeOme (keys %{$excludedOmes}) {
					$memberFail = 1 if ($member =~ m/^$excludeOme/);
				}
			}
			next if ($memberFail == 1);
			$member2group{$member} = $group;
		}
	}
	return(\%member2group);	
}

sub Parse_group_file_by_og {
	my ($clusteringOutfile, $excludedOmes) = @_;
	my $datestring = localtime();					
	my (%group2member);
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cant open $clusteringOutfile for reading, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		next if ($line =~ m/^Orthogroup/);
		my (@members) = split/\s+/, $line;
		my $group = shift @members;
		$group =~ s/:$//;
		foreach my $member (@members) {
			my $memberFail = 0;
			if (defined $excludedOmes) {
				foreach my $excludeOme (keys %{$excludedOmes}) {
					$memberFail = 1 if ($member =~ m/^$excludeOme/);
				}
			}
			next if ($memberFail == 1);
			push @{$group2member{$group}}, $member;
		}
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


sub Parse_repeat_GFF {
	my ($ome2gffFile) = @_;
	my %repeats;
	foreach my $ome (keys %{$ome2gffFile}) {
		open (my $IN, '<', $ome2gffFile->{$ome}) or die("Error: can't open $ome2gffFile->{$ome} for reading\n");
		while (my $line = <$IN>) {
			next if ($line =~ m/^#/);
			chomp $line;
			my ($contigID, $annotator, $featureType, $begin, $end) = split("\t", $line);
			next if (($end - $begin) < $MINREPEAT);
			$repeats{$contigID}{$begin} = 1;
			$repeats{$contigID}{$end} = 1;
		}
	}
	return(\%repeats);	
}

sub Parse_feature_coords {
	my ($bedFile) = @_;
	my %coords;
	open (my $IN, '<', $bedFile) or die("Error: can't open $bedFile for reading\n");
	while (my $line = <$IN>) {
		chomp $line;
		# DK001_Scaffold1	97800	100136	DK001_V000019	$tag	+	DK001_ship00001
		my ($contigID, $begin, $end, $featureID, $idtag, $strand, $regionIDString, $annotation) = split("\t", $line);
		push @{$coords{$featureID}{$contigID}}, $begin, $end;
	}
	return(\%coords);	
}


sub Opts_check {
	my ($opts) = @_;
	usage() if (defined $opts->{'h'});
	usage("\nError: please provide a file to -o\n") if (not defined $opts->{'o'});
	usage("\nError: the file provided to -o does not exist\n") if (! -f $opts->{'o'});
	usage("\nError: please provide a file to -b\n") if (not defined $opts->{'b'});
	usage("\nError: the file provided to -b does not exist\n") if (! -f $opts->{'b'});
	usage("\nError: the file provided to -r does not exist\n") if ((defined $opts->{'r'}) && (! -f $opts->{'r'}));
	usage("\nError: please provide an integer to -g\n") if (not defined $opts->{'g'});
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
	if (defined $opts->{'e'}) {
		usage("\nError: the file provided to -e does not exist\n") if (! -f $opts->{'e'});
	}
}
