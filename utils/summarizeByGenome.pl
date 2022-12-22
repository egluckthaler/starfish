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
usage: summarizeByGenome.pl [options]

Required:
-a, --assembly    FILE   2 column tsv: assembly name, path to assembly FASTA
-b, --bed         FILE   BED file with Starship feature coordinates (output by starfish summarize)
-f, --feat        FILE   starships.feat file (output by starfish summarize)
-s, --separator   STR    the character separating genomeID from featureID (default: '_')

Optional:
-g, --gff         FILE   2 column tsv: genome code, path to GFF
-n, --nameField   STR    GFF3 attribute field where gene features are named (default: 'Name=')
-h, --help              print more details and exit
/;
	if (not defined $message) {
		$message = qq/
For each familyID in --family, this script will summarize the number of starships per
genome that are in single copy, multi copy, or absent. Prints to STDOUT. \n\n/;
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
		'assembly|a=s',
		'feat|f=s',
		'separator|s=s',
		'gff|g=s',
		'nameField|namefield|n=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	# summarize genome size
	my ($assemblyPaths) = dim_1_hash($opts{'assembly'}, "\t", "0:1");
	my %ome2genomeSize;
	foreach my $omeID (keys %{$assemblyPaths}) {
		my ($contigs) = Fasta_hash($assemblyPaths->{$omeID});
		foreach my $contig (keys %{$contigs}) {
			$ome2genomeSize{$omeID}+= length($contigs->{$contig});
		}
	}
	
	# summarize total genes
	my (%ome2genomeGenes);
	if (defined $opts{'gff'}) {
	my ($gffPaths) = dim_1_hash($opts{'gff'}, "\t", "0:1");
		foreach my $omeID (keys %{$gffPaths}) {
			my %countedGenes;
			my ($gffInfo) = Gff_sortable_gene_hash($gffPaths->{$omeID}, $opts{'nameField'});
			foreach my $contig (keys %{$gffInfo}) {
				foreach my $pos (keys %{$gffInfo->{$contig}}) {
					foreach my $seqID (keys %{$gffInfo->{$contig}->{$pos}}) {
						next if (exists $countedGenes{$seqID});
						$ome2genomeGenes{$omeID}++;
						$countedGenes{$seqID} = 1;
					}
				}
			}
		}
	}
	
	# parse starship lengths and counts, ignoring any nested starships
	my (%ome2starshipLength, %ome2starshipCount);
	my ($starship2length) = dim_1_hash($opts{'feat'}, "\t", "1:5");
	my ($nestedStarshipsTemp) = dim_0_hash($opts{'feat'}, "\t", "20");
	my %nestedStarships;
	foreach my $starshipString (keys %{$nestedStarshipsTemp}) {
		my (@starships) = split/,/, $starshipString;
		foreach my $starship (@starships) {
			$nestedStarships{$starship} = 1;
		}	
	}
	foreach my $starshipID (keys %{$starship2length}) {
		next if (exists $nestedStarships{$starshipID});
		my ($omeID) = split/$opts{'separator'}/, $starshipID;
		$ome2starshipLength{$omeID} += $starship2length->{$starshipID};
		$ome2starshipCount{$omeID}++;
	}
	
	# parse starshipgenes
	my ($ome2starshipGenes) = Parse_starships($opts{'bed'}, $opts{'separator'});
	
	# print results
	print "genomeID\tstarshipCount\tgenomeLength\tstarshipLength\tpercStarshipLength\tgeneCount\tstarshipGeneCount\tpercStarshipGenes\n";
	foreach my $omeID (sort keys %{$assemblyPaths}) {
		print "$omeID";
		if (exists $ome2starshipCount{$omeID}) {
			print "\t$ome2starshipCount{$omeID}";
		} else {
			print "\t0";
		}
		print "\t$ome2genomeSize{$omeID}";
		if (exists $ome2starshipLength{$omeID}) {
			print "\t$ome2starshipLength{$omeID}\t";
			print sprintf("%.3f", $ome2starshipLength{$omeID} / $ome2genomeSize{$omeID});
		} else {
			print "\t0\t0";
		}
		if (exists $ome2genomeGenes{$omeID}) {
			print "\t$ome2genomeGenes{$omeID}";
			if (exists $ome2starshipGenes->{$omeID}) {
				print "\t$ome2starshipGenes->{$omeID}\t";
				print sprintf("%.3f", $ome2starshipGenes->{$omeID} / $ome2genomeGenes{$omeID} );
			} else {
				print "\t0\t0";
			}
		} else {
				print "\tna\tna\tna";
		}
		print "\n";
	}
}

sub Parse_starships {
	my ($bedFile, $SEP) = @_;
	my (%ome2starshipCount, %examinedGenes, %ome2starshipGenes);
	open (my $IN, '<', $bedFile) or usage("\n\nerror: can't open $bedFile for reading\n");
	while (my $line = <$IN>) {
		chomp $line;
		# DK001_Scaffold1	97800	100136	DK001_V000019	$tag	+	DK001_s00001	ann
		my ($contigID, $begin, $end, $featureID, $idtag, $strand, $regionIDString, $annString) = split("\t", $line);
		my ($omeID) = split/$SEP/, $contigID;
		my (@regionIDs) = split/,/, $regionIDString; # to accomodate overlapping regions e.g., in case of nested starships
		# make sure to only look at genes that havent been counted before
		next if (exists $examinedGenes{$featureID} || $idtag eq 'insert' || $idtag eq 'flank' || $idtag eq 'extend');
		$ome2starshipGenes{$omeID}++;
		$examinedGenes{$featureID} = 1;
	}
	return(\%ome2starshipGenes);	
}

sub dim_1_hash{
my $usage = qq/
usage: my (\%hash) = dim_1_hash(<abs_path_to_file>, <separator>, <ordered_colon_separated_column_numbers>)
sub dim_1_hash inputs the path to a tab delimited file, and returns a hash where key = first column and value = second column.\n\n/;
	my ($file, $sep, $columnstring) = @_;
	die($usage, "Error: incorrect number of args..\n") if (scalar @_ != 3);
	die($usage, "Error: column numbers must be colon delimited\n") if ($columnstring !~ m/:/);
	my @coi = split/:/, $columnstring; #columns of interest
	die($usage, "Error: only two column numbers should be specified\n") if (scalar @coi != 2);
	my %hash;
	open(my $in, '<', $file) or die($usage, "Error: cannot open $file..\n");
	while (my $line = <$in>) {
		if ($line =~ m/^#/) {
			next;
		} else {
			chomp $line;
			my @cols = split/$sep/, $line;
			my ($first, $second) = ($cols[$coi[0]], $cols[$coi[1]]);
			die($usage, "Error: cannot parse line from $file:\n$line\n") if (not defined $first);
			die($usage, "Error: cannot parse line from $file:\n$line\n") if (not defined $second);
			$hash{$first} = $second;
		}
	}
	return(\%hash);
}

sub dim_0_hash{
my $usage = qq/
usage: my (\%hash) = dim_1_hash(<abs_path_to_file>, <separator>, <ordered_colon_separated_column_numbers>)
sub dim_1_hash inputs the path to a tab delimited file, and returns a hash where key = first column and value = second column.\n\n/;
	my ($file, $sep, $columnstring) = @_;
	die($usage, "Error: incorrect number of args..\n") if (scalar @_ != 3);
	my @coi = split/:/, $columnstring; #columns of interest
	my %hash;
	open(my $in, '<', $file) or die($usage, "Error: cannot open $file..\n");
	while (my $line = <$in>) {
		if ($line =~ m/^#/) {
			next;
		} else {
			chomp $line;
			my @cols = split/$sep/, $line;
			my ($first) = ($cols[$coi[0]]);
			die($usage, "Error: cannot parse line from $file:\n$line\n") if (not defined $first);
			$hash{$first} = 1;
		}
	}
	return(\%hash);
}

sub dim_2_hash{
my $usage = qq/
usage: my (\%hash) = dim_2_hash(<abs_path_to_file>, <separator>, <ordered_colon_separated_column_numbers>)
sub dim_2_hash inputs the path to a tab delimited file, and returns a hash where first key = first column, second key = second column and value = 1.\n\n/;
	my ($file, $sep, $columnstring) = @_;
	die($usage, "Error: incorrect number of args..\n") if (scalar @_ != 3);
	die($usage, "Error: column numbers must be colon delimited\n") if ($columnstring !~ m/:/);
	my @coi = split/:/, $columnstring; #columns of interest
	die($usage, "Error: only two column numbers should be specified\n") if (scalar @coi != 2);
	my %hash;
	open(my $in, '<', $file) or die($usage, "Error: cannot open $file..\n");
	while (my $line = <$in>) {
		if ($line =~ m/^#/) {
			next;
		} else {
			chomp $line;
			my @cols = split/$sep/, $line;
			my ($first, $second) = ($cols[$coi[0]], $cols[$coi[1]]);
			die($usage, "Error: cannot parse line from $file:\n$line\n") if (not defined $first);
			die($usage, "Error: cannot parse line from $file:\n$line\n") if (not defined $second);
			$hash{$first}{$second} = 1;
		}
	}
	return(\%hash);
}

sub Gff_sortable_gene_hash {
	my ($gffFile, $NAMEFIELD) = @_;
	my (%allGenes);
	open(my $IN, '<', $gffFile) or die("\nError: could not open $gffFile for reading\n");
	while (my $line = <$IN>) {
		last if ($line =~ m/^##FASTA/);
		next if ($line =~ m/^#/);
		chomp $line;
		#acquire info
		my ($contigID, $annotator, $featureType, $begin, $end, $NULL1, $strand, $NULL2, $info) = split("\t", $line);
		next if ($featureType ne 'gene'); # only parse gene features
		$info = "\t$info";
		$info =~ m/[;\t]$NAMEFIELD([^;]+)/;
		my $seqID = $1;
		if (not defined $seqID) {
			my $datestring = localtime();
			warn("[$datestring] warning: could not parse sequenceID from $line using the provided namefield separator $NAMEFIELD\n");
		} else {
			$allGenes{$contigID}{$begin}{$seqID} = "$seqID\t$line"; # to conform to structure of %neighborhood hashes
		}
	}
	my $datestring = localtime();
	die("[$datestring] warning: could not parse any 'gene' features from $gffFile, but they need to exist for this program to continue\n") if (scalar keys %allGenes == 0);
	return(\%allGenes);
}


sub Fasta_hash {
	my ($fa) = @_;
	my (%fasta, %examined_seqs);
	my ($header, $faFH);
	die("Error: the fasta file $fa supplied to Fasta_hash does not exit, exiting..\n") if (! -f $fa);
	if ($fa =~ /.gz$/) {
		open($faFH, "gunzip -c $fa |") or warn("Error: could not open pipe to $fa with gunzip\n") && return;
	} else {
		open ($faFH, '<', $fa) or warn("Error: could not open $fa\n") && return;
	}

	while (my $line = <$faFH>) {
		chomp $line;
		next if ($line =~ m/^\s/);
		if ($line =~ m/^>/) {
			$examined_seqs{$header} = 1 unless (not defined $header); #store the previous sequence's header
			$line =~ m/^>(.+)$/;
			$header = $1;
		} elsif (eof $faFH) { #ensures that the last header sequence entry will be loaded into hash
			$fasta{$header} .= $line if (not exists $examined_seqs{$header}); #ensures that sequences with same header are not added twice
		} else {
			$fasta{$header} .= $line if (not exists $examined_seqs{$header}); #ensures that sequences with same header are not added twice
		} 
	}
	return(\%fasta);
}



sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a file to --assembly\n") if (not defined $opts->{'assembly'});
	usage("\nError: the file provided to --assembly does not exist\n") if (! -f $opts->{'assembly'});
	usage("\nError: please provide a file to --bed\n") if (not defined $opts->{'bed'});
	usage("\nError: the file provided to --bed does not exist\n") if (! -f $opts->{'bed'});
	usage("\nError: please provide a file to --feat\n") if (not defined $opts->{'feat'});
	usage("\nError: the file provided to --feat does not exist\n") if (! -f $opts->{'feat'});
	usage("\nError: please provide a file to --assembly\n") if (not defined $opts->{'assembly'});
	usage("\nError: the file provided to --assembly does not exist\n") if (! -f $opts->{'assembly'});
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
	if (defined $opts->{'gff'}) {
		usage("\nError: the file provided to --gff does not exist\n") if (! -f $opts->{'gff'});
	}
	if (not defined $opts->{'nameField'}) {
		$opts->{'nameField'} = 'Name=';
	}
}
