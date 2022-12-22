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
usage: mmseqs2mclFormat.pl [options]

Required:
-i, --input         FILE   two column tsv with: mmseqs group ID, geneID
-g, --groupid       STR    a short string used as a prefix to name groups
-o, --outdir        DIR    the output directory

Optional:
-m, --mapping       FILE   two column tsv with: oldGeneID, newGeneID
-h, --help                 print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script reformats the output of a mmseqs easy-clust .tsv file to a mcl-formatted
groups file. New group ids will be assigned on the fly with prefix --group. 

if a file is provided to --matching, geneIDs in --input will be renamed according to their
new ID in --mapping, and a file with suffix '.idmap.mcl' will be printed, in addition to
a file with all original geneIDs, for record keeping. Useful for renaming captainIDs
according to their starshipID\n\n/;
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
		'groupid|g=s',
		'mapping|m=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	my $datestring = localtime();
	print "\n[$datestring] reading in data..\n";

	# parse file prefix for outfile
	my ($filePrefix, $null, $suffix) = fileparse($opts{'input'}, (".txt", ".tsv"));
	my $outfile = "$opts{'outdir'}/$filePrefix.mcl";
	my $outmapfile = "$opts{'outdir'}/$filePrefix.idmap.mcl" if (defined $opts{'mapping'});

	# parse old2newfile file
	my ($old2new) = dim_1_hash($opts{'mapping'}, "\t", "0:1") if (defined $opts{'mapping'});

	# parse og2gene file
	my ($og2gene) = dim_2_hash($opts{'input'}, "\t", "0:1");

	# organize by number of members
	my %og2size;
	foreach my $og (keys %{$og2gene}) {
		my $groupSize = 0;
		foreach my $geneID (keys %{$og2gene->{$og}}) {
			$groupSize++;
		}
		$og2size{$og} = $groupSize;
	}
	
	# print
	my ($geneCounter,$ogCounter) = (0,0);
	my ($OUT) = Open_FH($outfile);
	my ($OUTMAP) = Open_FH($outmapfile) if (defined $opts{'mapping'});
	my %observedIDs;
	foreach my $og (sort {$og2size{$b} <=> $og2size{$a}} keys %og2size) { # notice we iterate through og2size
		$ogCounter++;
		my $newOGid = sprintf("%04d", $ogCounter);
		if ($og2size{$og} != 1) {
			print $OUT "${opts{groupid}}${newOGid}:";
			print $OUTMAP "${opts{groupid}}${newOGid}:" if (defined $old2new);
		} else {
			print $OUT "sng$newOGid:";
			print $OUTMAP "sng$newOGid:" if (defined $old2new);
		}
		foreach my $geneID (sort keys %{$og2gene->{$og}}) {
			
			# for unmapped id file
			print $OUT "\t$geneID";

			# for mapped id file
			if (defined $opts{'mapping'}) {
				my $newGeneID = $geneID;
				$newGeneID = $old2new->{$geneID} if (exists $old2new->{$geneID});
				if (not exists $observedIDs{$newGeneID}) {
					print $OUTMAP "\t$newGeneID" if (defined $old2new);
					$observedIDs{$newGeneID} = 1;
				} else {
					$datestring = localtime();
					print "[$datestring] warning: the new geneID $newGeneID has been observed >1 time, printing old gene ID $geneID to preserve mutual exclusivity of groups\n";
					print $OUTMAP "\t$geneID" if (defined $old2new);
				}
			}
			$geneCounter++;
		}
		print $OUT "\n";
		print $OUTMAP "\n" if (defined $old2new);
	}
	
	$datestring = localtime();
	print "[$datestring] sorted $geneCounter sequences into ".scalar(keys %{$og2gene})." ortholog groups..\n";

	$datestring = localtime();
	print "[$datestring] done\n";
	
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

sub dim_2_hash{
	my ($file, $sep, $columnstring) = @_;
	usage("\nError: incorrect number of args to dim_2_hash..\n") if (scalar @_ != 3);
	usage("\nError: column numbers must be colon delimited for dim_2_hash\n") if ($columnstring !~ m/:/);
	my @coi = split/:/, $columnstring; #columns of interest
	usage("\nError: only two column numbers should be specified for dim_2_hash\n") if (scalar @coi != 2);
	my %hash;
	open(my $in, '<', $file) or usage("\nError: cannot open $file..\n");
	while (my $line = <$in>) {
		if ($line =~ m/^#/) {
			next;
		} else {
			chomp $line;
			my @cols = split/$sep/, $line;
			my ($first, $second) = ($cols[$coi[0]], $cols[$coi[1]]);
			usage("\nError: cannot parse line from $file:\n$line\n") if (not defined $first);
			usage("\nError: cannot parse line from $file:\n$line\n") if (not defined $second);
			$hash{$first}{$second} = 1;
		}
	}
	return(\%hash);
}

sub dim_1_hash{
	my ($file, $sep, $columnstring) = @_;
	my @coi = split/:/, $columnstring; #columns of interest
	my %hash;
	open(my $in, '<', $file) or usage("Error: cant open $file..\n");
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

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a short prefix for naming groups to --groupid\n") if (not defined $opts->{'groupid'});
	usage("\nError: please provide an the path to the output directory\n") if (not defined $opts->{'outdir'});
	usage("\nError: the provided output directory does not exist\n") if (! -d $opts->{'outdir'});
	usage("\nError: please provide a mmseqs2gene file for --input\n") if (not defined $opts->{'input'});
	usage("\nError: the mmseqs2gene file provided to --input does not exist\n") if (! -f $opts->{'input'});
	if (defined $opts->{'mappings'}) {
		usage("\nError: the mapping file provided to --mappings does not exist\n") if (! -f $opts->{'mappings'});
	}
}