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
usage: geneOG2mclFormat.pl [options]

Required:
-i, --input         FILE   two column tsv with: geneID, ortholog group ID
-o, --outdir        DIR    the output directory

Optional:
-r, --reverse              indicates --input is a tsv with: ortholog group ID, geneID
-h, --help                 print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script reformats a gene2OG file to a mcl-formatted groups file\n\n/;
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
		'reverse|r',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	my $datestring = localtime();
	print "\n[$datestring] reading in data..\n";

	# parse file prefix for outfile
	my ($filePrefix, $null, $suffix) = fileparse($opts{'input'}, (".txt"));
	my $outfile = "$opts{'outdir'}/$filePrefix.mcl";

	# parse gene2og file
	my $og2gene;
	if (defined $opts{'reverse'}) {
		($og2gene) = dim_2_hash($opts{'input'}, "\t", "0:1");
	} else {
		($og2gene) = dim_2_hash($opts{'input'}, "\t", "1:0");	
	}
	
	# print
	my $geneCounter = 0;
	my ($OUT) = Open_FH($outfile);
	foreach my $og (sort keys %{$og2gene}) {
		print $OUT "$og:";
		foreach my $geneID (sort keys %{$og2gene->{$og}}) {
			print $OUT "\t$geneID";
			$geneCounter++;
		}
		print $OUT "\n";
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

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide an the path to the output directory\n") if (not defined $opts->{'outdir'});
	usage("\nError: the provided output directory does not exist\n") if (! -d $opts->{'outdir'});
	usage("\nError: please provide a gene2OG file for --input\n") if (not defined $opts->{'input'});
	usage("\nError: the gene2OG file provided to --input does not exist\n") if (! -f $opts->{'input'});
}