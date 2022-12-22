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
-i    FILE     input file (accepted formats: fasta, tsv)
-r    FILE     2 column tsv: search term, replacement\n/;

	if (not defined $message) {
		$message = qq/
This script replaces all search terms in the user provided files with the corresponding
replacement terms. All search terms must be whole words, and only fields ending in \t or :
will be replaced in -i. Replaced text will be printed to STDOUT. \n\n/; 
	} else {
		$message = "\n$message\nuse -h for more details\n" ;
	}	
	die($usage, $message);
}

main: {
	my %opts;
	getopt('r:i:h', \%opts);
	Opts_check(\%opts);
	my ($termFile, $inputFile) = ($opts{'r'}, $opts{'i'});
	
	#load up search and replacement terms by genome
	my ($term2replace) = dim_1_hash($termFile, "\t", "0:1");
	
	#replace all instances of search term with corresponding replacement term in file
	open (my $IN, '<', $inputFile) or usage("Error: cannae open $inputFile for reading\n");
	while (my $line = <$IN>) {
		chomp $line;
		my @fields = split/\t/, $line;
		my $maxIndex = (scalar @fields) - 1;
		my $currentIndex = -1;
		foreach my $field (@fields) {
			$currentIndex++;
			if ($field =~ m/^>/) {
				$field =~ s/>//; #replace for fasta files
				if (exists $term2replace->{$field}) {
					print ">$term2replace->{$field}";
				} else {
					print ">$field";
				}
			} else {
				
				# replace all instances of x with y
				foreach my $x (keys %{$term2replace}) {
					my $y = $term2replace->{$x};
					$field =~ s/$x([:]*)/$y$1/g;
				}
				
				# strict replacement rules
# 				if (exists $term2replace->{$field}) {
# 					print "$term2replace->{$field}";
# 				} else {
# 					print "$field";
# 				}
			}
			print "$field";
			print "\t" if ($currentIndex != $maxIndex);
		}
		print "\n";
	}
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
	usage() if (defined $opts->{'h'});
	usage("Error: please provide a 2 column tsv listing search and replacement terms to -r\n") if (not defined $opts->{'r'});
	usage("Error: can't find the 2 column tsv provided to -r\n") if (! -f $opts->{'r'});
	usage("Error: please provide a path to the input file to -i\n") if (not defined $opts->{'i'});
	usage("Error: can't find the input file provided to -i\n") if (! -f $opts->{'i'});
}