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
usage: filterID.pl [options]

Required:
-i, --input         FILE   .txt or .fasta file in which to edit sequence IDs
-s, --separator     STR    the character separating genomeID from featureID (default: '_')
-o, --outdir        DIR    the output directory

Optional:
-f, --fields        STR    comma-separated list of columns in --input to edit (only for .txt file; 0-indexed; incompatible with -n)
-n, --namefield     STR    the GFF3 attribute field where sequence features are named (only for .gff file; incompatible with -f)
-h, --help                 print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script is meant to patch a temporary bug in glofish consolidate where concatenated
sequence IDs sometimes have the same sequence ID specified multiple times. If no --fields
are specificied, attempts sequence ID editing on each line\n\n/;
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
		'input|i=s',
		'namefield|nameField|n=s',
		'fields|f=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	my $datestring = localtime();
	print "\n[$datestring] reading in data..\n";
	
	# parse fields
	my @fields;
	@fields = split/,/, $opts{'fields'} if (defined $opts{'fields'});
	
	# read and print out simultaneously
	my ($filename) = fileparse($opts{'input'});
	my $outfile = "$opts{'outdir'}/$filename.seqedit";
	open (my $IN, '<', $opts{'input'}) or usage("error: can't open $opts{input} for reading\n");
	open (my $OUT, '>', $outfile) or usage("error: can't open $outfile for writing\n");
	my $SEP = $opts{'separator'};
	my $SEPprint = $SEP;
	$SEPprint =~ s/\\//g;
	my $NAMEFIELD = $opts{'namefield'} if (defined $opts{'namefield'});

	while (my $line = <$IN>) {
		chomp $line;
		my (@columns) = split/\t/, $line;
		if (scalar @fields > 0) {
			foreach my $field (@fields) {
				if (defined $columns[$field]) {
					my ($featureID, $info1, $info2) = split/\|/, $columns[$field];
					my ($genomeCode, $seqID) = split/$SEP/, $featureID;
					if ((defined $genomeCode) && (defined $seqID)) {
						my (@parts) = split/:/, $seqID;
						my %newParts;
						foreach my $part (@parts) {
							$newParts{$part} = 1;
						}
						my $newID = "$genomeCode$SEPprint".join(":", sort keys %newParts);
						$newID = "$newID|$info1|$info2" if (defined $info1);
						$columns[$field] = $newID;
					}
				}
			}
		} elsif ((defined $NAMEFIELD) && (defined $SEP)) {
			my $attributes = $columns[8];
			$attributes = "\t$attributes";
			$attributes =~ m/[;\t]$NAMEFIELD([^;]+)/;
			my $featureID = $1;
			my ($genomeCode, $seqID) = split/$SEP/, $featureID;
			my $newID;
			if ((defined $genomeCode) && (defined $seqID)) {
				my (@parts) = split/:/, $seqID;
				my %newParts;
				foreach my $part (@parts) {
					$newParts{$part} = 1;
				}
				$newID = "$genomeCode$SEPprint".join(":", sort keys %newParts);
			}
			$columns[8] =~ s/$NAMEFIELD$featureID/$NAMEFIELD$newID/g;
		} else {
			my $index = 0;
			foreach my $column (@columns) {
				my ($featureID, $info1, $info2) = split/\|/, $column;
				my ($genomeCode, $seqID) = split/$SEP/, $featureID;
				if ((defined $genomeCode) && (defined $seqID)) {
					my (@parts) = split/:/, $seqID;
					my %newParts;
					foreach my $part (@parts) {
						$newParts{$part} = 1;
					}
					my $newID = "$genomeCode$SEPprint".join(":", sort keys %newParts);
					$newID = "$newID|$info1|$info2" if (defined $info1);
					$columns[$index] = $newID;
				}
				$index++;
			}
		}
		print $OUT join("\t", @columns)."\n";
	}

	$datestring = localtime();
	print "[$datestring] done\n";
	
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide an the path to the output directory\n") if (not defined $opts->{'outdir'});
	usage("\nError: the provided output directory does not exist\n") if (! -d $opts->{'outdir'});
	usage("\nError: please provide the path to an input file to --input\n") if (not defined $opts->{'input'});
	usage("\nError: the file provided to --input does not exist\n") if (! -f $opts->{'input'});
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