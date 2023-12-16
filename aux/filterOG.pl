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
usage: filterOG.pl [options]

Required:
-a, --absent        INT    the maximum number of genomes missing any given orthogroup
-c, --copies        INT    the maximum number of copies of any given orthogroup in any genome
-O, --orthologs     FILE   a MCL-formatted file of all ortholog groups (e.g., output by orthofinder)
-o, --outdir        DIR    the output directory

Required, with defaults:
-s, --separator     STR    the character separating genomeID from featureID (default: '_')

Optional:
-h, --help               print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script filters a MCL-formatted orthogroups file for orthogroups that meet user-
defined criteria based on the maximum number of genomes missing the orthogroup (--absent)
and the maximum number of copies of any given orthogroup in any genome (--copies).

A filtered file containining only ortholog groups that meet the above criteria will be
printed to a file with the same prefix as --orthologs to --outdir\n\n/;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($usage, $message);
}

main: {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'outdir|o=s',
		'separator|s=s',
		'orthologs|O=s',
		'absent|a=i',
		'copies|c=i',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	my $datestring = localtime();
	print "\n[$datestring] reading in data..\n";

	# parse file prefix for outfile
	my ($filePrefix, $null) = fileparse($opts{'orthologs'}, (".mcl", ".txt"));
	my $outfile = "$opts{'outdir'}/$filePrefix.a$opts{'absent'}.c$opts{'copies'}.txt";
	
	# parse orthogroups
	my ($group2members) = Parse_group_file_by_og($opts{'orthologs'});
	
	# parse all genome codes present in orthogroups file
	my ($genomeCodes) = Parse_genomes_with_orthogroups($group2members, $opts{'separator'});
	
	my $filteredOGs;

	$datestring = localtime();
	print "[$datestring] filtering ".scalar(keys %{$group2members})." OGs..\n";

	###########################
	#### FILTER FOR ABSENT ####
	###########################
	
	($filteredOGs) = Filter_by_absence($group2members, $genomeCodes, $opts{'separator'}, $opts{'absent'});
	
	###########################
	#### FILTER FOR COPIES ####
	###########################

	($filteredOGs) = Filter_by_copies($filteredOGs, $opts{'separator'}, $opts{'copies'});
	
	#######################
	#### PRINT RESULTS ####
	#######################
	
	$datestring = localtime();
	print "[$datestring] found ".scalar(keys %{$filteredOGs})." OGs absent from max $opts{absent} out of ".scalar(keys %{$genomeCodes})." genomes with max $opts{copies} copies in any given genome\n";

	Print_results($filteredOGs, $outfile);

	$datestring = localtime();
	print "[$datestring] done\n";
	
}

sub Print_results {
	my ($group2members, $outfile) = @_;
	my ($OUT) = Open_FH($outfile);
	
	foreach my $groupID (sort keys %{$group2members}) {
		print $OUT "$groupID: ";
		print $OUT join(" ", @{$group2members->{$groupID}});
		print $OUT "\n";
	}
}

sub Filter_by_copies {
	my ($group2members, $SEP, $MAXCOPIES) = @_;
	my %filteredOGs;
	
	foreach my $groupID (keys %{$group2members}) {
		my %genome2copy;
		foreach my $member (@{$group2members->{$groupID}}) {
			my ($omeCode, $featureID) = split/$SEP/, $member;
			usage("\nError: can't parse a genome code from the gene '$member' in the orthogroups file using the separator '$SEP'\n") if (not defined $omeCode || not defined $featureID);
			$genome2copy{$omeCode}++;
		}
		my $copyFail = 0;
		foreach my $genomeCode (keys %genome2copy) {
			$copyFail = 1 if ($genome2copy{$genomeCode} > $MAXCOPIES);
			last if ($copyFail == 1);
		}
		push @{$filteredOGs{$groupID}}, @{$group2members->{$groupID}} if ($copyFail == 0);
	}
	return(\%filteredOGs);	
}


sub Filter_by_absence {
	my ($group2members, $genomeCodes, $SEP, $MAXABSENT) = @_;
	my %filteredOGs;
	
	foreach my $groupID (keys %{$group2members}) {
		my %genomesPresent;
		foreach my $member (@{$group2members->{$groupID}}) {
			my ($omeCode, $featureID) = split/$SEP/, $member;
			usage("\nError: can't parse a genome code from the gene '$member' in the orthogroups file using the separator '$SEP'\n") if (not defined $omeCode || not defined $featureID);
			$genomesPresent{$omeCode} = 1;
		}
		push @{$filteredOGs{$groupID}}, @{$group2members->{$groupID}} if ((scalar(keys %{$genomeCodes}) - scalar(keys %genomesPresent)) <= $MAXABSENT);
	}
	return(\%filteredOGs);	
}

sub Parse_genomes_with_orthogroups {
	my ($group2members, $SEP) = @_;
	my %genomeCodes;
	foreach my $groupID (keys %{$group2members}) {
		foreach my $member (@{$group2members->{$groupID}}) {
			my ($omeCode, $featureID) = split/$SEP/, $member;
			usage("\nError: can't parse a genome code from the gene '$member' in the orthogroups file using the separator '$SEP'\n") if (not defined $omeCode || not defined $featureID);
			$genomeCodes{$omeCode} = 1;
		}
	}
	return(\%genomeCodes);	
}

sub Parse_group_file_by_og {
	my ($clusteringOutfile) = @_;
	my $datestring = localtime();					
	my %group2member;
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cannot read $clusteringOutfile, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		my (@members) = split/\s+/, $line;
		my $group = shift @members;
		$group =~ s/:$//;
		push @{$group2member{$group}}, @members;
	}
	return(\%group2member);	
}

sub Open_FH {
my $usage = qq/
my (<fh>) = Open_FH(<filename>, <header>)
sub Open_FH opens a file for writing and prints out a header, if supplied\n\n/;
	my ($file, $header) = @_;
	open(my $fileout, '>', $file) or die($usage, "Error: cannae open $file for writing\n");
	$fileout->autoflush(1);
	print $fileout "$header" if (defined $header);
	return($fileout);
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide an the path to the output directory\n") if (not defined $opts->{'outdir'});
	usage("\nError: the provided output directory does not exist\n") if (! -d $opts->{'outdir'});
	usage("\nError: please provide a MCL-formatted file listing all ortholog groups in your genomes of interest for --orthologs\n") if (not defined $opts->{'orthologs'});
	usage("\nError: the file listing all ortholog groups provided to --orthologs does not exist\n") if (! -f $opts->{'orthologs'});
	if (not defined $opts->{'absent'}) {
		usage("\nError: please specify --absent, the maximum number of genomes missing any given orthogroup\n");
	} elsif ($opts->{'absent'} !~ m/^\d+/) {
		usage("\nError: the --absent argument must be an integer\n");
	}
	if (not defined $opts->{'copies'}) {
		usage("\nError: please specify --copies, the maximum number of copies of any given orthogroup in any genome\n");
	} elsif ($opts->{'copies'} !~ m/^\d+/) {
		usage("\nError: the --copies argument must be an integer\n");
	}
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