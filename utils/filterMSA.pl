#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(:config auto_abbrev no_ignore_case);
use FileHandle;
$|=1;

# check dependencies
my @commandlines = ("hhsearch");
Commandline_check(\@commandlines);

# default hhsearch parameters
my $hhsearchDefaults = "-M first -Z 250 -z 1 -B 0 -b 0 -v 1";

sub usage {
	my $message = shift;
	my $usage = qq/
usage: filterMSA.pl [options]

Required:
-i, --input         FILE   multiple sequence alignment in FASTA format
-m, --matches       STR    comma-separated list of database identifiers for filtering sequences
-o, --outdir        DIR    the output directory

Required, with defaults:
-d, --database      FILE   path to a hhsuite database (default: \/data1\/emile\/pdb70\/pdb70)
-a, --args          STR    quoted string of additional arguments to pass to hhsearch (default: '')

Optional:
-r, --restrict      FILE   single column of sequence IDs to restrict analysis to
-h, --help                 print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script executes hhsearch for each query sequence in a multiple sequence alignment. A
modified MSA containing only sequences with hits to user-specified database matches
(--matches) will then be printed to --outdir as --input.hhfilt.

WARNING: many files are temporarily printed to --outdir\n\n/;
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
		'args|a=s',
		'database|d=s',
		'matches|m=s',
		'restrict|r=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	my $datestring = localtime();
	print "\n[$datestring] reading in data..\n";

	# parse file prefix for outfile
	my ($filePrefix, $null, $suffix) = fileparse($opts{'input'}, (".mafft", ".fas", ".fasta", ".faa", ".fna", ".clipkit"));
	my $outfile = "$opts{'outdir'}/$filePrefix.hhfilt$suffix";

	# parse fasta file, and keep position of sequence in file
	my ($seqs) = Fasta_hash_conserve_pos($opts{'input'});
	
	# parse filterseqs, if available
	my ($restrictIDs) = dim_0_hash($opts{'restrict'}, "\t", "0") if (defined $opts{'restrict'}); 
	
	# parse matches
	my @matches = split/,/, $opts{'matches'};
	
	$datestring = localtime();
	if (not defined $restrictIDs) {
		print "[$datestring] filtering ".scalar(keys %{$seqs})." sequences in $opts{input}..\n";
	} else {
		print "[$datestring] filtering ".scalar(keys %{$restrictIDs})." sequences in $opts{input}..\n";
	}
	##########################
	#### EXECUTE HHSEARCH ####
	##########################
	
	my %hhrFiles;
	
	foreach my $pos (sort {$a <=> $b} keys %{$seqs}) {
		my ($queryOUT) = Open_FH("$opts{outdir}/query.fasta");

		# there will only ever be 1 sequence per pos
		my ($headerSpace, $headerNospace);
		foreach my $header (keys %{$seqs->{$pos}}) {
			$headerSpace = $header;
			$header =~ m/^([^\s]+)/;
			$headerNospace = $1;
			print $queryOUT ">$header\n$seqs->{$pos}->{$header}\n";
		}
		
		if ((not defined $restrictIDs) || ((exists $restrictIDs->{$headerSpace}) || (exists $restrictIDs->{$headerNospace}))) {
		
			# now print out all other sequence in alignment, below the first focal query sequence
			foreach my $pos2 (sort {$a <=> $b} keys %{$seqs}) {
				# skip the focal sequence, obvs
				next if ($pos eq $pos2);
			
				# there will only ever be 1 sequence per pos
				foreach my $header2 (keys %{$seqs->{$pos2}}) {
					print $queryOUT ">$header2\n$seqs->{$pos2}->{$header2}\n";
				}
			}
		
			my $hhrOutfile = "$opts{outdir}/$headerNospace.hhr";
			$hhrFiles{$hhrOutfile} = $headerSpace;
		
			# skip hhsearch execution if output file exists
			if (! -f $hhrOutfile) {
				my ($failCheck) = system("hhsearch -i $opts{outdir}/query.fasta -d $opts{database} -o $hhrOutfile $hhsearchDefaults $opts{args}");
				$datestring = localtime();					
				if ($failCheck != 0) { warn "\n\n[$datestring] error: could not execute hhsearch on commandline for $headerNospace, exiting..\n$!\n";}
			} else {
				$datestring = localtime();
				print "[$datestring] $hhrOutfile already exists, skipping hhsearch of $headerNospace..\n";
			}
		}
		system("rm $opts{outdir}/query.fasta");
	}
	
	###############################
	#### FILTER SEARCH RESULTS ####
	###############################

	my %filteredSeqs;
	
	foreach my $hhrFile (keys %hhrFiles) {
		my $matchFound = 0;
		$datestring = localtime();
		open(my $IN, '<', $hhrFile) or usage("\n\n[$datestring] error: cannot read $hhrFile, exiting\n");
		while (my $line = <$IN>) {
			chomp $line;
			$line =~ s/ +/\t/g;
			my ($null, $no, $hitID) = split/\t/, $line;
			next if (not defined $hitID);
			foreach my $match (@matches) {
				$matchFound = 1 if ($hitID eq $match);
			}
		}
		my $ogHeader = $hhrFiles{$hhrFile};
		$filteredSeqs{$ogHeader} = 1 if ($matchFound == 1);
	}


	#######################
	#### PRINT RESULTS ####
	#######################
	
	$datestring = localtime();
	print "[$datestring] kept ".scalar(keys %filteredSeqs)." sequences with hits to --matches\n";
	my $filtOUT = Open_FH("$outfile");

	foreach my $pos (sort {$a <=> $b} keys %{$seqs}) {
		foreach my $header (keys %{$seqs->{$pos}}) {
			print $filtOUT ">$header\n$seqs->{$pos}->{$header}\n" if (exists $filteredSeqs{$header});
		}
	}
		
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

sub Fasta_hash_conserve_pos {
	my ($fa) = @_;
	my (%fasta, %examined_seqs);
	my ($header, $faFH);
	die("Error: the fasta file $fa supplied to Fasta_hash does not exit, exiting..\n") if (! -f $fa);
	if ($fa =~ /.gz$/) {
		open($faFH, "gunzip -c $fa |") or warn("Error: could not open pipe to $fa with gunzip\n") && return;
	} else {
		open ($faFH, '<', $fa) or warn("Error: could not open $fa\n") && return;
	}
	my $counter = 0;
	while (my $line = <$faFH>) {
		chomp $line;
		next if ($line =~ m/^\s/);
		if ($line =~ m/^>/) {
			$examined_seqs{$header} = 1 unless (not defined $header); #store the previous sequence's header
			$line =~ m/^>(.+)$/;
			$header = $1;
			$counter++;
		} elsif (eof $faFH) { #ensures that the last header sequence entry will be loaded into hash
			$fasta{$counter}{$header} .= $line if (not exists $examined_seqs{$header}); #ensures that sequences with same header are not added twice
		} else {
			$fasta{$counter}{$header} .= $line if (not exists $examined_seqs{$header}); #ensures that sequences with same header are not added twice
		} 
	}
	return(\%fasta);
}

sub Commandline_check {
	my ($args) = @_;
	my $datestring = localtime();
	foreach my $arg (@{$args}) {
# 		my ($failcheck) = `$arg 2> /dev/null`;
# 		($failcheck) = `$arg -h 2> /dev/null` if (not defined $failcheck);
# 		($failcheck) = `$arg -help 2> /dev/null` if (not defined $failcheck);
# 		($failcheck) = `$arg --help 2> /dev/null` if (not defined $failcheck);
		my $failcheck = `which $arg 2>&1`;
		if (($failcheck =~ m/no $arg/) || ($failcheck =~ m/^$/)) { die("[$datestring] error: can't exec $arg on commandline (are you sure its in your PATH?), exiting..\n");}
	}
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
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a comma-separated list of database accessions to --matches\n") if (not defined $opts->{'matches'});
	usage("\nError: please provide an the path to the output directory\n") if (not defined $opts->{'outdir'});
	usage("\nError: the provided output directory does not exist\n") if (! -d $opts->{'outdir'});
	usage("\nError: please provide a multiple sequence alignment in FASTA format for --input\n") if (not defined $opts->{'input'});
	usage("\nError: the multiple sequence alignment provided to --orthologs does not exist\n") if (! -f $opts->{'input'});
	usage("\nError: the file provided to --restrict does not exist\n") if ((defined $opts->{'restrict'}) && (! -f $opts->{'restrict'}));
	usage("\nError: please provide at least 1 database identifier to --matches for filtering sequences\n") if (not defined $opts->{'matches'});
	if (not defined $opts->{'database'}) {
		print "\nno hhsuite database defined, defaulting to \/data1\/emile\/pdb70\/pdb70\n";
		$opts->{'database'} = "/data1/emile/pdb70/pdb70";
	}
	if (not defined $opts->{'args'}) {
		print "\nno additional arguments for hhsearch specified, defaulting to none";
		$opts->{'args'} = "";
	}
}