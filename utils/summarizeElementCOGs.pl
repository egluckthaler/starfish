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
usage: summarizeElementCOGs.pl [options]

Required:
-q, --queries     FILE   three column tsv: geneID, starshipID, functionalCategory
-r, --report      FILE   blastp output in 'outfmt 6' format
-b, --bed         FILE   BED file with Starship features
-o, --omes        FILE   one column tsv: omeID
-s, --separator   STR    the character separating genomeID from featureID (default: '_')

Optional:
-h, --help              print more details and exit
/;
	if (not defined $message) {
		$message = qq/
For each genome in --omes and each gene in --queries, this script will determine whether
the gene's best hit (by bit score) is in a starship or the genomic background or missing,
and then will summarize the total counts for each functional category. Genes that are best
hits to multiple query genes will only be counted once to contribute to a particular
functional category count. Currently customized for COG categories, includes function to
parse combined COG annotations. Prints to STDOUT. \n\n/;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($usage, $message);
}

my %cog2info = (
	"D" => "Cellular processes and signaling\tCell cycle control, cell division, chromosome partitioning",
	"M" => "Cellular processes and signaling\tCell wall/membrane/envelope biogenesis",
	"N" => "Cellular processes and signaling\tCell motility",
	"O" => "Cellular processes and signaling\tPost-translational modification, protein turnover, and chaperones",
	"T" => "Cellular processes and signaling\tSignal transduction mechanisms",
	"U" => "Cellular processes and signaling\tIntracellular trafficking, secretion, and vesicular transport",
	"V" => "Cellular processes and signaling\tDefense mechanisms",
	"W" => "Cellular processes and signaling\tExtracellular structures",
	"Y" => "Cellular processes and signaling\tNuclear structure",
	"Z" => "Cellular processes and signaling\tCytoskeleton",
	"A" => "Information storage and processing\tRNA processing and modification",
	"B" => "Information storage and processing\tChromatin structure and dynamics",
	"J" => "Information storage and processing\tTranslation, ribosomal structure and biogenesis",
	"K" => "Information storage and processing\tTranscription",
	"L" => "Information storage and processing\tReplication, recombination and repair",
	"C" => "Metabolism\tEnergy production and conversion",
	"E" => "Metabolism\tAmino acid transport and metabolism",
	"F" => "Metabolism\tNucleotide transport and metabolism",
	"G" => "Metabolism\tCarbohydrate transport and metabolism",
	"H" => "Metabolism\tCoenzyme transport and metabolism",
	"I" => "Metabolism\tLipid transport and metabolism",
	"P" => "Metabolism\tInorganic ion transport and metabolism",
	"Q" => "Metabolism\tSecondary metabolism",
	"R" => "Unknown\tGeneral function prediction only",
	"S" => "Unknown\tFunction unknown");

main: {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'bed|b=s',
		'queries|q=s',
		'report|r=s',
		'omes|o=s',
		'separator|s=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	# parse query info
	my ($queryGene2starship) = dim_1_hash($opts{'queries'}, "\t", "0:1");
	my ($starship2queryGene) = dim_2_hash($opts{'queries'}, "\t", "1:0");
	my ($queryGene2category) = dim_1_hash($opts{'queries'}, "\t", "0:2");
	my ($starshipIDs) = dim_0_hash($opts{'queries'}, "\t", "1");
	
	# parse omes
	my ($omeIDs) = dim_0_hash($opts{'omes'}, "\t", "0");
	
	# parse all genes in starships, and consider any gene is a starship with prefix _s as being in a full starship
	# and any gene in a starship with prefix _e as being in a partial starship
	my ($gene2starship) = Parse_genes_in_starships($opts{'bed'});

	# parse best hit per gene query per genome in the blast output, according to bitscore
	# structured: {genomeID}{functionalCategory}{compartmentType}++
	my ($query2counts) = Parse_hits_per_genome($queryGene2starship, $gene2starship, $queryGene2category, $omeIDs, $opts{'report'}, $opts{'separator'});

	# print results
	print "genomeID\tCOGcategory\tCOGdescription\thitsInBoundedStarships\thitsInPartialStarships\thitsInAllStarships\thitsInBackground\thitsMissing\n";
	foreach my $omeID (sort keys %{$query2counts}) {
		foreach my $category (sort keys %{$query2counts->{$omeID}}) {
			print "$omeID\t$cog2info{$category}";
			my $totalPresent = 0;
			if (exists $query2counts->{$omeID}->{$category}->{'starships'}) {
				print "\t$query2counts->{$omeID}->{$category}->{'starships'}";
				$totalPresent+=$query2counts->{$omeID}->{$category}->{'starships'};
			} else {
				print "\t0";
			}
			if (exists $query2counts->{$omeID}->{$category}->{'partial'}) {
				print "\t$query2counts->{$omeID}->{$category}->{'partial'}";
				$totalPresent+=$query2counts->{$omeID}->{$category}->{'partial'};
			} else {
				print "\t0";
			}
			print "\t$totalPresent";
			if (exists $query2counts->{$omeID}->{$category}->{'background'}) {
				print "\t$query2counts->{$omeID}->{$category}->{'background'}";
				$totalPresent+=$query2counts->{$omeID}->{$category}->{'background'};
			} else {
				print "\t0";
			}
			if (exists $query2counts->{$omeID}->{$category}->{'missing'}) {
				print "\t$query2counts->{$omeID}->{$category}->{'missing'}";
				$totalPresent+=$query2counts->{$omeID}->{$category}->{'missing'};
			} else {
				print "\t0";
			}
			print "\n";
		}
	}
}

sub Parse_hits_per_genome {
	my ($queryGene2starship, $gene2starship, $queryGene2category, $omeIDs, $blastReport, $SEP) = @_;
	open(my $in, '<', $blastReport) or usage("Error: cannae open $blastReport for reading\n");
	
	# find the best hit per genome
	# structured: {genomeID}{functionalCategory}{compartmentType}++
	my (%query2hits, %genomeMax);
	while (my $line = <$in>) {
		chomp $line;
		my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split/\t/, $line;
		if (exists $queryGene2starship->{$qseqid}) { # make sure query gene is in a query Starship
			my ($omeCode) = split/$SEP/, $sseqid;
			if (not exists $genomeMax{$omeCode}{$qseqid}) {
				$genomeMax{$omeCode}{$qseqid} = $bitscore;
				$query2hits{$omeCode}{$qseqid} = $sseqid;
			} elsif ($bitscore > $genomeMax{$omeCode}{$qseqid}) {
				$genomeMax{$omeCode}{$qseqid} = $bitscore;
				$query2hits{$omeCode}{$qseqid} = $sseqid;
			}
		}
	}
	
	# count the number of genomes with either _s, _e or neither
	my (%query2cat, %observedBestHits);
	foreach my $queryID (keys %{$queryGene2starship}) { # iterate through each query gene to check if it has a hit or not in each genome
		foreach my $omeID (keys %{$omeIDs}) {
			my @categories = split//, $queryGene2category->{$queryID}; # for cog categories that are sometimes combined
			foreach my $category (@categories) {
				if (exists $query2hits{$omeID}{$queryID}) { # was a best hit recovered?
					my $bestHit = $query2hits{$omeID}{$queryID};
					next if (exists $observedBestHits{$bestHit}); # avoid double counting hits in genome that hit to multiple query genes
					$observedBestHits{$bestHit} = 1;
					if (exists $gene2starship->{$bestHit}) { # is the hit in a starship?
						if ($gene2starship->{$bestHit} =~ m/_s/) {
							$query2cat{$omeID}{$category}{'starships'}++;
						} elsif ($gene2starship->{$bestHit} =~ m/_e/) {
							$query2cat{$omeID}{$category}{'partial'}++;
						}
					} else {
						$query2cat{$omeID}{$category}{'background'}++;
					}
				} else {
					$query2cat{$omeID}{$category}{'missing'}++;
				}
			}
		}
	}
	return(\%query2cat);	
}

sub Parse_genes_in_starships {
	my ($bedFile) = @_;
	my %gene2starship;
	open (my $IN, '<', $bedFile) or usage("\n\nerror: can't open $bedFile for reading\n");
	while (my $line = <$IN>) {
		chomp $line;
		# DK001_Scaffold1	97800	100136	DK001_V000019	$tag	+	DK001_s00001	ann
		my ($contigID, $begin, $end, $featureID, $idtag, $strand, $regionIDString, $annString) = split("\t", $line);
		my (@regionIDs) = split/,/, $regionIDString; # to accomodate overlapping regions e.g., in case of nested starships
		foreach my $regionID (@regionIDs) {
			$gene2starship{$featureID} = $regionID;
		}
	}
	return(\%gene2starship);	
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

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a file to --bed\n") if (not defined $opts->{'bed'});
	usage("\nError: the file provided to --bed does not exist\n") if (! -f $opts->{'bed'});
	usage("\nError: please provide a file to --queries\n") if (not defined $opts->{'queries'});
	usage("\nError: the file provided to --queries does not exist\n") if (! -f $opts->{'queries'});
	usage("\nError: please provide a file to --report\n") if (not defined $opts->{'report'});
	usage("\nError: the file provided to --report does not exist\n") if (! -f $opts->{'report'});
	usage("\nError: please provide a file to --omes\n") if (not defined $opts->{'omes'});
	usage("\nError: the file provided to --omes does not exist\n") if (! -f $opts->{'omes'});
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


