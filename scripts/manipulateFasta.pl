#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
$| = 1;

sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]

Required:
-f    STR     fastafile input
-o    STR     output directory

Optional:
-n    STR     name of output file (not applicable to -d or -p; will override default output prefix)
-p    INT     maximum number of sequences per output file
-d    STR     delimeter for parsing genome code from sequence ID to print to genome-specific file (incompatible with -p)
-l    INT     minimum sequence length (in residues; smaller will be excluded)
-r    STR     residues to exclude (e.g., -)
-x    STR     path to file listing names of sequences to exclude
-i    STR     path to file listing names of sequences to include

manipulate_fasta.pl is a general purpose tool to manipulate multifasta files. Users can
break up a fasta file into smaller files containing fewer sequences; can exlude sequences
based on a minimum sequence length; can exclude or include sequences based on their name;
can break the fasta up by genome code, as specified by a delimiter (first delimited field
will be selected as genome code)

Headers in output file are modified to contain only first string of non-whitespace characters\n\n/;
	die($usage, $message);
}
 
main: {

	my %opts;
	getopt('f:p:o:l:x:i:d:r:n:h', \%opts);
	Opts_check(\%opts);
	my ($fastafile, $maxseq, $minlength, $delimiter, $excludefile, $includefile, $restrictedResidue, $outdir) = ($opts{'f'}, $opts{'p'}, $opts{'l'}, $opts{'d'}, $opts{'x'}, $opts{'i'}, $opts{'r'},$opts{'o'});
	my ($filename, $fastadir, $suffix) = fileparse($fastafile, (".fasta", ".fa", ".fna", ".faa"));
	my ($fasta) = Fasta_hash_nospace($fastafile);
	
	#filter for size
	($fasta) = Size_filter($fasta, $minlength) if (defined $minlength);
	
	#filter by excluded sequences
	($fasta) = Name_filter_exclude($fasta, $excludefile) if (defined $excludefile);
	
	#filter by included sequences
	($fasta) = Name_filter_include($fasta, $includefile) if (defined $includefile);

	#filter out restricted residues
	($fasta) = Residue_filter($fasta, $restrictedResidue) if (defined $restrictedResidue);

	#for breaking up fasta by max seq
	if (defined $maxseq) {
		my ($part, $seqcount) = (0, 0);
		my $OUT;
		foreach my $header (sort keys %{$fasta}) {
			if ($seqcount == 0 || $seqcount >= $maxseq) {
				$part++;
				$seqcount = 1; #reset seqcount, or initialize it if seq = 0
				($OUT) = Open_FH("$outdir/$filename.filt$part$suffix");
				$header =~ m/^([^\s]+)/; #grab first string of nonwhitespace characters
				my $modheader = $1;
				if (exists $fasta->{$header}) {
					$fasta->{$header} =~ s/\*//g;
					print $OUT ">$modheader\n$fasta->{$header}\n";
				} else {
					warn("Warning: can't find sequence associated with $header, skipping printout\n");
				}
			} else {
				$header =~ m/^([^\s]+)/; #grab first string of nonwhitespace characters
				my $modheader = $1;
				$fasta->{$header} =~ s/\*//g;
				if (exists $fasta->{$header}) {
					print $OUT ">$modheader\n$fasta->{$header}\n";
					$seqcount++;
				} else {
					warn("Warning: can't find sequence associated with $header, skipping printout\n");
				}
			}
		}
		print "$part fasta files have been created with at most $maxseq sequences each\n";

	#for breaking up fasta by genome
	} elsif (defined $delimiter) {
		my $genomeCount = 0;
		my ($genomeFasta) = Genome_filter($fasta, $delimiter);
		foreach my $genomeCode (sort keys %{$genomeFasta}) {
			$genomeCount++;
			my $outfile = "$outdir/$genomeCode$suffix";
			my ($OUT) = Open_FH($outfile);
			foreach my $header (sort keys %{$genomeFasta->{$genomeCode}}) {
				$header =~ m/^([^\s]+)/; #grab first string of nonwhitespace characters
				my $modheader = $1;
				$genomeFasta->{$genomeCode}->{$header} =~ s/\*//g;
				print $OUT ">$modheader\n$genomeFasta->{$genomeCode}->{$header}\n";
			}
		}
		print "$genomeCount fasta files have been created with genome-specific sequences\n";
	} else { #otherwise print out normally
		my $OUT;
		if (defined $opts{'n'}) {
			($OUT) = Open_FH("$outdir/$opts{'n'}"); 
		} else {
			($OUT) = Open_FH("$outdir/$filename.filt$suffix");
		}
		foreach my $header (sort keys %{$fasta}) {
			$header =~ m/^([^\s]+)/; #grab first string of nonwhitespace characters
			my $modheader = $1;
			$fasta->{$header} =~ s/\*//g;
			print $OUT ">$modheader\n$fasta->{$header}\n";
		}
	}
}

sub Genome_filter {
	my ($fasta, $SEP) = @_;
	my %genomeFasta;
	foreach my $header (keys %{$fasta}) {
		my ($genomeCode) = split/$SEP/, $header;
		if (defined $genomeCode) { 
			$genomeFasta{$genomeCode}{$header} = $fasta->{$header};
		} else {
			warn("Warning: can't parse a genome code using delimiter $SEP for $header, skipping\n");
		}
	}
	return(\%genomeFasta);	
}

sub Name_filter_exclude {
	my ($fasta, $excludefile) = @_;
	my %editfasta;
	my $excludecount = 0;
	my ($exclude) = dim_0_hash($excludefile, "\t", "0");
	foreach my $header (keys %{$fasta}) {
		if (exists $exclude->{$header}) {
			$excludecount++;
		} else {
			$editfasta{$header} = $fasta->{$header};
		}
	}
	print "$excludecount out of ".scalar(keys %{$exclude})." sequences listed in $excludefile have been removed\n";
	return(\%editfasta);	
}

sub Name_filter_include {
	my ($fasta, $includefile) = @_;
	my %editfasta;
	my $includecount = 0;
	my ($include) = dim_0_hash($includefile, "\t", "0");
	foreach my $header (keys %{$fasta}) {
		if (exists $include->{$header}) {
			$includecount++;
			$editfasta{$header} = $fasta->{$header};
		} 
	}
	print "$includecount out of ".scalar(keys %{$include})." sequences listed in $includefile have been retained\n";
	return(\%editfasta);	
}

sub Residue_filter {
	my ($fasta, $restrictedResidue) = @_;
	my %editfasta;
	my $excludecount = 0;
	foreach my $header (keys %{$fasta}) {
		$fasta->{$header} =~ s/$restrictedResidue//g;
		$editfasta{$header} = $fasta->{$header};
	}
	return(\%editfasta);	
}

sub Size_filter {
	my ($fasta, $minlength) = @_;
	my %editfasta;
	my $excludecount = 0;
	foreach my $header (keys %{$fasta}) {
		if (length($fasta->{$header}) >= $minlength) {
			$editfasta{$header} = $fasta->{$header};
		} else {
			$excludecount++;
		}
	}
	print "$excludecount sequences < ${minlength}bp have been removed\n";
	return(\%editfasta);	
}

sub Open_FH {
	my ($file, $header) = @_;
	open(my $fileout, '>', $file) or usage("Error: can't open output file $file for writing\n");
	$fileout->autoflush(1);
	print $fileout "$header" if (defined $header);
	return($fileout);
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

sub Fasta_hash_nospace {
	my ($fa) = @_;
	usage("Error: the fasta file $fa does not exist\n") if (! -f $fa);
	my %fasta;
	my ($header, $faFH);
	if ($fa =~ /.gz$/) {
		open($faFH, "gunzip -c $fa |") or warn("Error: could not open pipe to $fa with gunzip\n") && return;
	} else {
		open ($faFH, '<', $fa) or warn("Error: could not open $fa\n") && return;
	}
	while (my $line = <$faFH>) {
		chomp $line;
		next if ($line =~ m/^\s/);
		if ($line =~ m/^>/) {
			$line =~ m/^>([^\s]+)/;
			$header = $1;
		} elsif (eof $faFH) { #ensures that the last header sequence entry will be loaded into hash
			$fasta{$header} .= $line;
		} else {
			$fasta{$header} .= $line;
		} 
	}
	return(\%fasta);
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (defined $opts->{'h'});
	usage("Error: please provide a fasta file to -f\n") if (not defined $opts->{'f'});
	usage("Error: can't find the fasta file provided to -f\n") if (! -f $opts->{'f'});
	usage("Error: please specify an output directory to -o\n") if (not defined $opts->{'o'});
	usage("Error: can't find the output directory provided to -o\n") if (! -d $opts->{'o'});
	if (defined $opts->{'p'}) {
		usage("Error: the argument for -p must be an integer\n") if ($opts->{'p'} !~ m/^\d+$/);
	}
	if (defined $opts->{'l'}) {
		usage("Error: the argument for -l must be an integer\n") if ($opts->{'l'} !~ m/^\d+$/);
	}
	if (defined $opts->{'x'}) {
		usage("Error: can't find the file provided to -x\n") if (! -f $opts->{'x'});
	}
	if (defined $opts->{'i'}) {
		usage("Error: can't find the file provided to -i\n") if (! -f $opts->{'i'});
	}
	if (defined $opts->{'d'}) {
		if ($opts->{'d'} eq ':') {
			usage("\nError: the delimiter character provided to -d cannot be ':'\n");
		} elsif ($opts->{'d'} eq ';') {
			usage("\nError: the delimiter character provided to -d cannot be ';'\n");
		} elsif ($opts->{'d'} eq '|') {
			usage("\nError: the delimiter character provided to -d cannot be '|'\n");
		}
		if (defined $opts->{'d'}) {
			$opts->{'d'} = quotemeta($opts->{'d'});	# to allow splitting on special characters, like '.'
		}
		if ((defined $opts->{'p'})) {
			usage("\nError: you cannot specify both -d and -p");
		}
	}
}
