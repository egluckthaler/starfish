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
usage: summarizeFamilyByGenome.pl [options]

Required:
-f, --family      FILE   two column tsv: starshipID, familyID
-b, --bed         FILE   BED file with Starship features
-o, --omes        FILE   two column tsv: omeID, category of interest
-s, --separator   STR    the character separating genomeID from featureID (default: '_')

Optional:
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
		'family|f=s',
		'omes|o=s',
		'separator|s=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################

	# parse query info
	my ($starship2family) = dim_1_hash($opts{'family'}, "\t", "0:1");
	my ($familyIDs) = dim_0_hash($opts{'family'}, "\t", "1");
	
	# parse omes
	my ($omeIDs) = dim_1_hash($opts{'omes'}, "\t", "0:1");
	
	# parse starships
	# {omeID}{starshipID} = 1
	my ($genome2starship) = Parse_starships($opts{'bed'},$opts{'separator'});
	
	# sum family per genome
	my ($genome2family) = Sum_across_family($starship2family, $genome2starship);
	
	# my
	print "genomeID\tcategory\tfamilyID\tcountCategory\tvalue\n";
	foreach my $omeID (sort keys %{$omeIDs}) {
		foreach my $familyID (sort keys  %{$familyIDs}) {
			print "$omeID\t$omeIDs->{$omeID}\t$familyID\t";
			if (exists $genome2family->{$omeID}->{$familyID}) {
				if ($genome2family->{$omeID}->{$familyID} == 1) {
					print "singleCopy\t1\n";
				} else {
					print "multiCopy\t1\n";
				}
			} else {
				print "absent\t1\n";		
			}
		}
	}
}

sub Sum_across_family {
	my ($starship2family, $genome2starship) = @_;
	my %ome2family;
	foreach my $omeID (keys %{$genome2starship}) {
		foreach my $starshipID (keys %{$genome2starship->{$omeID}}) {
			if (exists $starship2family->{$starshipID}) { # only look at starships that are assigned to a family of interest
				$ome2family{$omeID}{$starship2family->{$starshipID}}++;
			}
		}
	}
	return(\%ome2family);	
}

sub Parse_starships {
	my ($bedFile, $SEP) = @_;
	my %ome2starship;
	open (my $IN, '<', $bedFile) or usage("\n\nerror: can't open $bedFile for reading\n");
	while (my $line = <$IN>) {
		chomp $line;
		# DK001_Scaffold1	97800	100136	DK001_V000019	$tag	+	DK001_s00001	ann
		my ($contigID, $begin, $end, $featureID, $idtag, $strand, $regionIDString, $annString) = split("\t", $line);
		my ($omeID) = split/$SEP/, $contigID;
		my (@regionIDs) = split/,/, $regionIDString; # to accomodate overlapping regions e.g., in case of nested starships
		foreach my $regionID (@regionIDs) {
			$ome2starship{$omeID}{$regionID} = 1;
		}
	}
	return(\%ome2starship);	
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
	usage("\nError: please provide a file to --family\n") if (not defined $opts->{'family'});
	usage("\nError: the file provided to --family does not exist\n") if (! -f $opts->{'family'});
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





