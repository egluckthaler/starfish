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
usage: og2ann.pl [options]

Required:
-a, --ann     FILE   two column tsv with: geneID, annotation
-g, --group   FILE   a MCL-formatted orthogroups file

Optional:
-h, --help           print more details and exit
/;
	if (not defined $message) {
		$message = qq/
This script is for annotating orthogroups. It reports the frequency of each annotation 
assigned to the set of sequences in each OG. Prints to STDOUT. \n\n/;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($usage, $message);
}

main : {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'ann|a=s',
		'group|g=s',
		'h|help');
	Opts_check(\%opts);

	######################
	#### READ IN DATA ####
	######################
	
	# Retrieve seqs associated with OGoi
	my ($prot2og) = Parse_orthogroups_file($opts{'group'});

	# Retrieve annotations for seqs
	my ($prot2ann) = Parse_annotation_file($opts{'ann'});
			
	########################
	#### PARSE ANN FREQ ####
	########################

	# Assess annotation freq
	my ($OGannfreq) = Parse_annotation_freq($prot2ann, $prot2og);

	#######################
	#### PRINT RESULTS ####
	#######################

	foreach my $og (sort {$OGannfreq->{$a} <=> $OGannfreq->{$b}} keys %{$OGannfreq}) {
		print "$og\t";
		foreach my $ann (sort { $OGannfreq->{$og}->{$b} <=> $OGannfreq->{$og}->{$a} } keys %{$OGannfreq->{$og}}) {
			print "$ann\t$OGannfreq->{$og}->{$ann});";
		}
		print "\n";
	}
}

sub Parse_annotation_file {
	my ($annotationfile) = @_;
	my %prot2ann;
	open(my $annin, '<', $annotationfile) or usage("Error: can't open $annotationfile for reading\n");
	while(my $line = <$annin>) {
		next if ($line=~ m/^#/);
		chomp $line;
		my ($prot, $annotation) = split/\t/, $line;
		push @{$prot2ann{$prot}}, $annotation;
	}
	return(\%prot2ann);
}

sub Parse_annotation_freq {
	my ($prot2ann, $prot2og) = @_;
	my (%OGfreq, %protfreq);
	foreach my $prot (keys %{$prot2ann}) {
		if (exists $prot2og->{$prot}) { # make sure this prot is assigned to an OG
			$protfreq{$prot2og->{$prot}}{$prot} = 1;
			foreach my $ann (@{$prot2ann->{$prot}}) {
				#Use this parser if you don't care about splitting up comma separated descriptions
				$ann = lc $ann; #some annotations are annoyingly not controlled by case (e.g., you can have chitin and Chitin)
				$ann =~ s/^\s+//;
				$OGfreq{$prot2og->{$prot}}{$ann}++;
			}
		}
	}
	my (%annfreq);
	foreach my $og (keys %OGfreq) {
		my $protcount = scalar keys %{$protfreq{$og}}; #how many prots in this OG?
		foreach my $ann (keys %{$OGfreq{$og}}) {
			my $anncount = $OGfreq{$og}{$ann};
			my $freq = sprintf("%.3f", $anncount/$protcount); #divide the number of annotations by the number of prots in OG
			$annfreq{$og}{$ann} = $freq;
		}
	}
	return(\%annfreq);
}

sub Parse_orthogroups_file {
	my ($MCLFILE) = @_;
	my %seqs;
	open(my $mclin, '<', $MCLFILE) or usage("Error: can't open $MCLFILE for reading\n");
	while (my $line = <$mclin>) {
		chomp $line;
		my @cols = split/\s+/, $line;
		my $OG = shift @cols;
		$OG =~ s/:$//;
		my @seqs;
		foreach my $seq (@cols) {
			$seqs{$seq} = $OG;
		}
	}
	return(\%seqs);	
}

sub dim_0_hash{
my $usage = qq/
usage: my (\%hash) =  dim_0_hash(<abs_path_to_file>, <separator>, <column_numbers>)
sub 0_dim_hash inputs the path to a tab delimited file, and returns a hash where key = column number and value = 1\n\n/;
	my ($file, $sep, $coi) = @_;
	die($usage, "Error: incorrect number of args..\n") if (scalar @_ != 3);
	die($usage, "Error: only one column number should be specified\n") if ($coi =~ m/[^\d]/);
	my %hash;
	open(my $in, '<', $file) or die($usage, "Error: cannot open $file..\n");
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
	usage("\nError: no arguments provided\n") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a file to --ann\n") if (not defined $opts->{'ann'});
	usage("\nError: the file provided to --ann does not exist\n") if (! -f $opts->{'ann'});
	usage("\nError: please provide a file to --group\n") if (not defined $opts->{'group'});
	usage("\nError: the file provided to --group does not exist\n") if (! -f $opts->{'group'});
}

