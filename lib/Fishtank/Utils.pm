package Fishtank::Utils;
use strict;
use warnings;
use FileHandle;
use File::Basename;
use Data::Dumper;
use Text::Levenshtein qw(distance);
use Exporter qw(import);
our @EXPORT_OK = qw(
Add_genome_to_sequence
Reverse_gene_hash
Parse_upstream_region_from_gene
Filter_out_value
Filter_in_first_key
Filter_out_first_key
Parse_duf3435_from_regions
Reverse_complement
Fasta_hash
Fasta_hash_nospace
Fasta_hash_many_files
Gff_gene_hash
Gff_sortable_gene_hash
Glofish_bed_hash
Sort_bed_hash_by_coordinate
Annotation_hash
dim_0_hash
dim_1_hash
dim_2_hash
Open_FH
Commandline_check
Format_check_fasta
Format_check_fasta_paths
Format_check_gff
Format_name
Calculate_pairwise_mashdist
Calculate_mashdist_candidate_known
Homopolymer_check
Edit_distance_check
Parse_known_elements
Parse_region_boundaries

);

####### USEFUL UNIVERSAL SUBROUTINES ########

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

sub Homopolymer_check {
	# fail check if sequence is a homopolymer
	my ($seq) = @_;
	my $failCheck = 0;
	$seq = lc($seq);
	$failCheck = 1 if ($seq =~ m/^[a]+$|^[t]+$|^[g]+$|^[c]+$/);
	return($failCheck);		
}

sub Edit_distance_check {
	# fail check if sequences are less than minSim similar
	# minsim should range from (0-1]
	
	my ($seq1, $seq2, $minSim) = @_;
	die("Error: min similarity for sub Edit_distance check must range from (0,1]\n") if (($minSim <= 0) || ($minSim > 1));
	my $failCheck = 0;
	
	# Levenshtein edit distance
	my $lev = distance($seq1, $seq2);
	
	my $longestLength = length($seq1);
	$longestLength = length($seq2) if (length($seq2) > $longestLength);

	# convert to % similarity: p = 1 - l/m, Where l is the levenshtein distance and m is the length of the longest of the two words
	my $seqSim = 1 - ($lev / $longestLength);
	$failCheck = 1 if ($seqSim < $minSim);
	return($failCheck, $seqSim);
}

sub Add_genome_to_sequence {
	my ($seqs, $SEP) = @_;
	my %genomeSeqs;
	foreach my $header (keys %{$seqs}) {
		my ($genome) = split/$SEP/, $header;
		$genomeSeqs{$genome}{$header} = $seqs->{$header};
	}
	return(\%genomeSeqs);
}

sub Reverse_gene_hash {
	# input  : {contigID}{starshipID}{geneID} = [begin, end, strand, tag, annotation]
	# returns: {geneID}{starshipID}{contigID} = [begin, end, strand, tag, annotation]
	my ($genes) = @_;
	my %reversed;
	foreach my $contigID (keys %{$genes}) {
		foreach my $regionID (keys %{$genes->{$contigID}}) {
			foreach my $geneID (keys %{$genes->{$contigID}->{$regionID}}) {
				my ($begin, $end, $strand, $tag, $annotation) = @{$genes->{$contigID}->{$regionID}->{$geneID}};
				push @{$reversed{$geneID}{$regionID}{$contigID}}, $begin, $end, $strand, $tag;
				push @{$reversed{$geneID}{$regionID}{$contigID}}, $annotation if (defined $annotation);
			}
		}
	}
	return(\%reversed);
}

sub Parse_upstream_region_from_gene {
	# upstream regions have the same ID as their associated gene
	# returns: {contigID}{starshipID}{geneID} = [begin, end, strand, 'up']
	# assumes following structures:
	# genes: {contigID}{starshipID}{geneID} = [begin, end, strand, tag]
	# assemblies: {contigID} = seq
	# upstream range: min-max
	
	my ($genes, $assemblies, $upstreamRange) = @_;
	my %upstreams;
	my ($minOffset, $maxOffset) = split/-/, $upstreamRange;
	#print Dumper(keys%{$genes});die;
	foreach my $contigID (keys %{$genes}) {
		foreach my $regionID (keys %{$genes->{$contigID}}) {
			my ($upBegin, $upEnd);
			foreach my $geneID (keys %{$genes->{$contigID}->{$regionID}}) {
				my ($begin, $end, $strand, $tag) = @{$genes->{$contigID}->{$regionID}->{$geneID}};
				# ensure that upstream coordinates are actually on the contig
				# orient the upstream coordinates relative to the gene's orientation
				if ($strand eq '+') {
					$upBegin = $begin - $maxOffset;
					$upEnd = $begin - $minOffset;
				} elsif ($strand eq '-') {
					$upBegin = $end + $minOffset;
					$upEnd = $end + $maxOffset;
				} else {
					die("Error: can't parse orientation of gene=$geneID from contig=$contigID region=$regionID\n");
				}

				# make sure coordinates are actually on the contig
				# covers all possible conditions, whether cap is + or -
				if (defined $assemblies->{$contigID}) {
					my $contigLength = length($assemblies->{$contigID});
					$upBegin = 1 if ($upBegin < 0);
					$upEnd = $contigLength if ($upEnd > $contigLength);
					print "No coordinates found for: $contigID\t$regionID\t$geneID\n" if (not defined $contigLength);
				
					# upstream region does not exist if its end is located off of the contig, or if it starts at position 1
					# upstream region does not exist if its beginning is located off of the contig
					# as a sanity check, make sure the coordinates are correctly oriented, and are at least 10bp long
					my $minUpstreamLength = 10;
					($upBegin, $upEnd) = ($upEnd, $upBegin) if ($upBegin > $upEnd);
					my $upLength = $upEnd - $upBegin + 1;
					
					if (($upLength < $minUpstreamLength) || ($upEnd < 2) || ($upBegin > $contigLength)) {
						push @{$upstreams{$contigID}{$regionID}{$geneID}}, 'NA', 'NA', 'NA', 'up';
					} else { # upstream region is located on contig!	 
						push @{$upstreams{$contigID}{$regionID}{$geneID}}, $upBegin, $upEnd, $strand, 'up';
					}
				} else {
					my $datestring = localtime();
					print "[$datestring] warning: can't find sequence for $contigID in BED file, skipping..\n";
				}
			}
		}
	}
	return(\%upstreams);
}

sub Reverse_complement {	
	my ($seq, $id) = @_;
	$seq = reverse($seq);
	$seq = lc($seq);
	my (@bases) = split//, $seq;
	my $newseq;
	my $warning = 0;
	foreach my $base (@bases) {
		my $newbase;
		if ($base eq 'a') {
			$newbase = 't';
		} elsif ($base eq 't') {
			$newbase = 'a';
		} elsif ($base eq 'g') {
			$newbase = 'c';
		} elsif ($base eq 'c') {
			$newbase = 'g';
		} elsif ($base eq 'n') {
			$newbase = 'n';
		} elsif ($base eq '-') {
			$newbase = '-';
		}
		if (defined $newbase) {
			$newseq .= $newbase;
		} else {
			my $datestring = localtime();
			print "[$datestring] warning: could not reverse complement $base in $id, inserting n. Silencing warnings for this sequence\n" if ($warning == 0);
			$warning = 1;
			$newseq .= 'n';
		}
	}
	return($newseq);	
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

sub Fasta_hash_nospace {
	my ($fa) = @_;
	die("Error: the fasta file $fa supplied to Fasta_hash does not exit, exiting..\n") if (! -f $fa);
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

# no spaces are loaded into header
sub Fasta_hash_many_files { 
	my ($files) = @_; # experts a hash where each key is the path to a file
	my %fasta;
	foreach my $fa (keys %{$files}) {
		die("Error: the fasta file $fa supplied to Fasta_hash does not exit, exiting..\n") if (! -f $fa);
		my ($header, $faFH);
		if ($fa =~ /.gz$/) {
			open($faFH, "gunzip -c $fa |") or warn("Error: could not open pipe to $fa with gunzip\n") && return;
		} else {
			open ($faFH, '<', $fa) or warn("Error: could not open $fa\n") && return;
		}
		while (my $line = <$faFH>) {
			chomp $line;
			if ($line =~ m/^>/) {
				$line =~ m/^>([^\s]+)/;
				$header = $1;
			} elsif (eof $faFH) { #ensures that the last header sequence entry will be loaded into hash
				$fasta{$header} .= $line;
			} else {
				$fasta{$header} .= $line;
			} 
		}
	}
	return(\%fasta);
}

sub Gff_gene_hash {
	my ($gffFile, $NAMEFIELD) = @_;
	my %genes;
	die("Error: the GFF3 attribute field where gene features are named is not being passed correctly to sub GFF_hash\n") if (not defined $NAMEFIELD);
	open (my $IN, '<', $gffFile) or die("Error: can't open $gffFile for reading\n");
	my $datestring = localtime();
	print "[$datestring] parsing $gffFile and skipping all features not labeled 'gene'..\n";
	while (my $line = <$IN>) {
		chomp $line;
		my ($contigID, $annotator, $featureType, $begin, $end, $NULL1, $strand, $NULL2, $attributes) = split("\t", $line);
		next if ($featureType ne 'gene'); # only parse gene features
		$attributes = "\t$attributes";
		$attributes =~ m/[;\t]$NAMEFIELD([^;]+)/;
		my $seqID = $1;
		$genes{$contigID}{$seqID} = $line;
	}
	return(\%genes);	
}

sub Glofish_bed_hash {
	my ($bedFile) = @_;
	my %regions;
	open (my $IN, '<', $bedFile) or die("Error: can't open $bedFile for reading\n");
	while (my $line = <$IN>) {
		chomp $line;
		# DK001_Scaffold1	97800	100136	DK001_V000019	$tag	+	DK001_ship00001
		my ($contigID, $begin, $end, $featureID, $idtag, $strand, $regionIDString, $annotation) = split("\t", $line);
		my (@regionIDs) = split/,/, $regionIDString; # to accomodate overlapping regions e.g., in case of nested starships
		foreach my $regionID (@regionIDs) {
			push @{$regions{$contigID}{$regionID}{$featureID}}, $begin, $end, $strand, $idtag;
			push @{$regions{$contigID}{$regionID}{$featureID}}, $annotation if (defined $annotation);
		}
	}
	return(\%regions);	
}

sub Sort_bed_hash_by_coordinate {
	# intended for sorting output of Glofish_bed_hash
	my ($regions) = @_;
	my %sorted;
	foreach my $contigID (keys %{$regions}) {
		foreach my $regionID (keys %{$regions->{$contigID}}) {
			foreach my $featureID (keys %{$regions->{$contigID}->{$regionID}}) {
				my ($begin, $end, $strand, $idtag, $annotation) = @{$regions->{$contigID}->{$regionID}->{$featureID}};
				@{$sorted{$contigID}{$regionID}{$begin}{$featureID}} = @{$regions->{$contigID}->{$regionID}->{$featureID}};
			}
		}
	}
	return(\%sorted);
}

sub Annotation_hash {
	my ($annotationFile) = @_;
	my %annotations;
	my %observedAnnotations;
	open (my $IN, '<', $annotationFile) or die("Error: can't open $annotationFile for reading\n");
	while (my $line = <$IN>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($geneID, $field, $annotation) = split("\t", $line);
		if (not exists $observedAnnotations{$geneID}{$annotation}) {
			push @{$annotations{$geneID}{$field}}, $annotation;
			$observedAnnotations{$geneID}{$annotation} = 1;
		} 
	}
	return(\%annotations);	
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

sub Format_check_fasta {
	# note that the first argument can either be the output of Fasta_hash or a path to a fasta file
	my ($fastaFile, $SEP) = @_;
	die("Error: the separator character is not being passed correctly to sub Format_check_fasta\n") if (not defined $SEP);
	if (ref $fastaFile eq ref {}) { # check if provided variable is a hash
		foreach my $header (keys %{$fastaFile}) {
			my $datestring = localtime();
			die("[$datestring] error: sequenceID $header in the provided Fasta hash contains a '|' character, but should not for parsing purposes. Please remove all '|' from all sequence IDs\n") if ($header =~ m/\|/);
			my @components = split/$SEP/, $header;
			die("[$datestring] error: $header in the provided Fasta hash is being parsed into >2 components using separator \'$SEP\'. Make sure ALL sequence headers are formatted like <genomeID><separator><featureID>\n") if (scalar @components > 2);
			die("[$datestring] error: $header in the provided Fasta hash is being parsed into <2 components using separator \'$SEP\'. Make sure ALL sequence headers are formatted like <genomeID><separator><featureID>\n") if (scalar @components < 2);
		}
	} else {
		my ($seqs) = Fasta_hash_nospace($fastaFile);
		foreach my $header (keys %{$seqs}) {
			my $datestring = localtime();
			die("[$datestring] error: sequenceID $header in $fastaFile contains a '|' character, but should not for parsing purposes. Please remove all '|' from all sequence IDs\n") if ($header =~ m/\|/);
			my @components = split/$SEP/, $header;
			die("[$datestring] error: $header in $fastaFile is being parsed into >2 components using separator \'$SEP\'. Make sure ALL sequence headers are formatted like <genomeID><separator><featureID>\n") if (scalar @components > 2);
			die("[$datestring] error: $header in $fastaFile is being parsed into <2 components using separator \'$SEP\'. Make sure ALL sequence headers are formatted like <genomeID><separator><featureID>\n") if (scalar @components < 2);
		}
	}
}

sub Format_check_gff {
	my ($gffFile, $SEP, $NAMEFIELD) = @_;
	open (my $IN, '<', $gffFile) or die("Error: can't open $gffFile for reading\n");
	die("Error: the separator character $SEP is not defined in sub Format_check_gff\n") if (not defined $SEP);
	die("Error: the namefield string $NAMEFIELD is not defined in sub Format_check_gff\n") if (not defined $NAMEFIELD);
	
	while (my $line = <$IN>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($contigID, $annotator, $featureType, $begin, $end, $NULL1, $strand, $NULL2, $info) = split("\t", $line);
		my $datestring = localtime();
		die("[$datestring] error: $gffFile should be tab-separated, but likely is not\n") if (not defined $featureType);
		die("[$datestring] error: $gffFile should only contain 'gene' features but contains $featureType\n") if ($featureType ne 'gene'); # only parse gene features
		$info = "\t$info"; # for helping with regex parsing
		$info =~ m/[;\t]$NAMEFIELD([^;]+)/;
		my $header = $1;
		if (defined $header) {
			my @components = split/$SEP/, $header;
			die("[$datestring] error: $header on line $line in $gffFile is being parsed into >2 components using separator \'$SEP\'. Make sure ALL gene feature names are formatted like <genomeID><separator><featureID>\n") if (scalar @components > 2);
			die("[$datestring] error: $header on line $line in $gffFile is being parsed into <2 components using separator \'$SEP\'. Make sure ALL gene feature names are formatted like <genomeID><separator><featureID>\n") if (scalar @components < 2);
		} else {
			die("[$datestring] error: $line in $gffFile does not have a parse-able featureID using namefield \'$NAMEFIELD\'. Make sure ALL gene feature names are are stored in the attributes column like <namefield><geneName>\n");
		}
	}
}

sub Format_name {
	my ($name, $genomeID, $SEP) = @_;
	my $newName;
	my $SEPprint = $SEP;
	$SEPprint =~ s/\\//g;

	# check if separator is already in name
	if ($name =~ m/$SEP/) {
		# check if genome ID is recoverable in first SEP field, because who in their right mind would put it anywhere else?
		my @components = split/$SEP/, $name;
		if ($components[0] eq $genomeID) {
			# paste all remaining fields together, removing all SEPs if they exist
			# will reconstitute any names that are already formatted correctly
			my $firstField = shift @components;
			if (scalar @components < 1) { # in the rare case where there is no parseable feature ID after genomeID, skip it
				$newName = 'NA';
			} else {
				my $featureID = join("", @components);
				$newName = "${genomeID}${SEPprint}${featureID}";
			}
		} else { 
			#format the whole thing, after removing all SEPs
			my $featureID = join("", @components);
			$newName = "${genomeID}${SEPprint}${featureID}";
		}
	} else {
		# format the whole thing
		$newName = "${genomeID}${SEPprint}${name}";
	}
	
	# remove any and all ':', ';' and '|' from name
	$newName =~ s/:|;|\|//g;
	return($newName);
}

sub Parse_duf3435_from_regions {
	# will not parse any sequences with 'cap' tag
	my ($starships, $duf3435TAG) = @_;
 	my (%orderedTaggedGenes, %duf3435s);

 	# first identify all genes with duf3435 Tags, and store them according to their 5' coordinate
 	# also identify all candidate sequences tagged with 'cap': if a region already has a candidate captain, then we should not return any other tyrs form that region for downstream analysis
 	my %shipsWithCap;
 	foreach my $contigID (keys %{$starships}) {
 		#next unless ($contigID eq 'DK189_Scaffold10'); # for debugging
		foreach my $starshipID (keys %{$starships->{$contigID}}) {
 			foreach my $geneID (keys %{$starships->{$contigID}->{$starshipID}}) {
 				my ($currentBegin, $currentEnd, $strand, $tag, $annotation) = @{$starships->{$contigID}->{$starshipID}->{$geneID}};
				push @{$orderedTaggedGenes{$contigID}{$starshipID}{$currentBegin}}, $geneID, @{$starships->{$contigID}->{$starshipID}->{$geneID}} if ($tag eq $duf3435TAG);
				$shipsWithCap{$starshipID} = 1 if ($tag eq 'cap');
			}
		}
	} 

	# Each region can have up to two putative captains:
	#     in the 5'-3' direction, the most upstream DUF3435 sequence that is also on the + strand (if the most upstream DUF3435 is on the - strand, no putative captain is returned)
	#     in the 3'-5' direction, the most downstream DUF3435 sequence that is also on the - strand (if the most downstream DUF3435 is on the + strand, no putative captain is returned)
	
 	foreach my $contigID (keys %orderedTaggedGenes) {
 		foreach my $starshipID (keys %{$orderedTaggedGenes{$contigID}}) {
 			
 			# skip starships that already have a candidate captain
 			next if (exists $shipsWithCap{$starshipID});
 		
			# sorting by hash key is a kind of awkward way of identifying the most up- and down-stream genes, but whatever
			# 5'-3'
			foreach my $firstPosition (sort {$a <=> $b} keys %{$orderedTaggedGenes{$contigID}{$starshipID}}) {
 				my ($geneID, $currentBegin, $currentEnd, $strand, $tag, $annotation) = @{$orderedTaggedGenes{$contigID}{$starshipID}{$firstPosition}};
				if ($strand eq '+') {
					push @{$duf3435s{$contigID}{$starshipID}{$geneID}}, $currentBegin, $currentEnd, $strand, $tag; # only consider a DUF3435 to be a captain if it is the first DUF3435 encountered in the 5'-3' direction AND oriented correctly
					push @{$duf3435s{$contigID}{$starshipID}{$geneID}}, $annotation if (defined $annotation);
				}
				last;
			}
			# 3'-5'
			foreach my $lastPosition (sort {$b <=> $a} keys %{$orderedTaggedGenes{$contigID}{$starshipID}}) { # only consider a DUF3435 to be a captain if it is the first DUF3435 encountered in the 3'-5' direction AND oriented correctly
 				my ($geneID, $currentBegin, $currentEnd, $strand, $tag, $annotation) = @{$orderedTaggedGenes{$contigID}{$starshipID}{$lastPosition}};
				 if ($strand eq '-') {
					push @{$duf3435s{$contigID}{$starshipID}{$geneID}}, $currentBegin, $currentEnd, $strand, $tag;
					push @{$duf3435s{$contigID}{$starshipID}{$geneID}}, $annotation if (defined $annotation);
				}
				last;
			}
		}
	}
	return(\%duf3435s);	
}

sub Filter_in_first_key {
	my ($idsToTest, $idsToKeep) = @_;
	my %filtered;
	foreach my $id (keys %{$idsToTest}) {
		$filtered{$id} = $idsToTest->{$id} if (exists $idsToKeep->{$id});
	}
	return(\%filtered);	
}

sub Filter_out_first_key {
	my ($idsToTest, $idsToRemove) = @_;
	my %filtered;
	foreach my $id (keys %{$idsToTest}) {
		$filtered{$id} = $idsToTest->{$id} if (not exists $idsToRemove->{$id});
	}
	return(\%filtered);	
}

sub Filter_out_value {
	my ($idsToTest, $idsToRemove) = @_;
	my %filtered;
	foreach my $id (keys %{$idsToTest}) {
		$filtered{$id} = $idsToTest->{$id} if (not exists $idsToRemove->{$idsToTest->{$id}});
	}
	return(\%filtered);	
}

sub Calculate_pairwise_mashdist {
	my ($genes, $assemblies, $OUTDIR) = @_;
	my %distances;
	
	# only consider queries from different starship regions from the ref
	foreach my $refGeneID (keys %{$genes}) {
		foreach my $refRegionID (keys %{$genes->{$refGeneID}}) {
			foreach my $refContigID (keys %{$genes->{$refGeneID}->{$refRegionID}}) {
				my ($refBegin, $refEnd, $refStrand, $refTag) = @{$genes->{$refGeneID}->{$refRegionID}->{$refContigID}};
				next if ($refBegin eq 'NA'); # skip regions not actually located on contigs
				
				# name files
				my $refFile = "$OUTDIR/mashRef.fa";
				my $queFile = "$OUTDIR/mashQue.fa";
				
				# print out ref sequence
				my ($refOUT) = Open_FH($refFile);
				my ($queOUT) = Open_FH($queFile);
				
				# substr is 0 indexed, and length includes starting position
				my $refSeq = substr($assemblies->{$refContigID}, $refBegin - 1, $refEnd - $refBegin + 1);
				print $refOUT ">$refGeneID\n$refSeq\n";
				
				# iterate through all other genes from other starship regions and print to queFile
				foreach my $queGeneID (keys %{$genes}) {
					foreach my $queRegionID (keys %{$genes->{$queGeneID}}) {
						next if ($refRegionID eq $queRegionID);
						foreach my $queContigID (keys %{$genes->{$queGeneID}->{$queRegionID}}) {
							my ($queBegin, $queEnd, $queStrand, $queTag) = @{$genes->{$queGeneID}->{$queRegionID}->{$queContigID}};
							next if ($queBegin eq 'NA'); # skip regions not actually located on contigs
							my $queSeq = substr($assemblies->{$queContigID}, $queBegin - 1, $queEnd - $queBegin + 1);
							print $queOUT ">$queGeneID\n$queSeq\n";
						}
					}
				}
				
				# execute mash
				my $mashString = `mash dist -i $refFile $queFile 2> $OUTDIR/mash.err`;
				usage("Error: could not execute mash on the commandline for reference $refGeneID\n") if (not defined $mashString);
				my @mashLines = split/\n/, $mashString;
				foreach my $line (@mashLines) {
					chomp $line;
					next if ($line =~ m/^Sketching/);
					my ($ref, $que, $dist, $pvalue, $sharedHashes) = split/\t/, $line;
					$distances{$ref}{$que} = $dist;
				}
				system("rm $refFile $queFile $OUTDIR/mash.err");
			}
		}
	}
	return(\%distances);
}

sub Calculate_mashdist_candidate_known {
	my ($candidateCaptains, $knownElements, $bedFeatures, $assemblies, $capmash, $OUTDIR) = @_;
	my %distances;
	
	foreach my $refGeneID (keys %{$candidateCaptains}) {

		#next if ($refGeneID ne 'mp070_08458');

		foreach my $refRegionID (keys %{$candidateCaptains->{$refGeneID}}) {
			foreach my $refContigID (keys %{$candidateCaptains->{$refGeneID}->{$refRegionID}}) {
				my ($refBegin, $refEnd, $refStrand, $refTag) = @{$candidateCaptains->{$refGeneID}->{$refRegionID}->{$refContigID}};
				next if ($refBegin eq 'NA'); # skip regions not actually located on contigs
				
				# name files
				my $refFile = "$OUTDIR/mashRef.fa";
				my $queFile = "$OUTDIR/mashQue.fa";
				
				# print out ref sequence
				my ($refOUT) = Open_FH($refFile);
				my ($queOUT) = Open_FH($queFile);
				
				# substr is 0 indexed, and length includes starting position
				my $refSeq = substr($assemblies->{$refContigID}, $refBegin - 1, $refEnd - $refBegin + 1);
				print $refOUT ">$refGeneID\n$refSeq\n";
				
				# iterate through all known captainsand print to queFile
				foreach my $queGeneID (keys %{$knownElements}) {
					foreach my $queRegionID (keys %{$knownElements->{$queGeneID}}) {
						foreach my $queContigID (keys %{$knownElements->{$queGeneID}->{$queRegionID}}) {
							my ($queBegin, $queEnd, $queStrand, $queTag) = @{$bedFeatures->{$queContigID}->{$queRegionID}->{$queGeneID}}; # notice we are recovering the captain coordinates here; knownElements contain element coordinates
							next if ($queBegin eq 'NA'); # skip regions not actually located on contigs
							my $queSeq = substr($assemblies->{$queContigID}, $queBegin - 1, $queEnd - $queBegin + 1);
							print $queOUT ">$queGeneID\n$queSeq\n";
						}
					}
				}
				
				# execute mash
				my $mashString = `mash dist -i $refFile $queFile 2> $OUTDIR/mash.err`;
				usage("Error: could not execute mash on the commandline for reference $refGeneID\n") if (not defined $mashString);
				my @mashLines = split/\n/, $mashString;
				foreach my $line (@mashLines) {
					chomp $line;
					next if ($line =~ m/^Sketching/);
					my ($ref, $que, $dist, $pvalue, $sharedHashes) = split/\t/, $line;
					
					# skip loading up any comparisons > capmash
					next if ($dist > $capmash);
					
					$distances{$ref}{$que} = $dist;
				}
				system("rm $refFile $queFile $OUTDIR/mash.err");
			}
		}
	}
	return(\%distances);
}

sub Parse_known_elements {
	# identifies known elements based on 'cap' idtag associated with captain sequence
	# input  : {contigID}{regionID}{geneID} = [featureBegin, featureEnd, strand, tag]
	my ($genes) = @_;
	my %elements;
	my %uniqueElements;
	my $elementCount = 0;
	# identify the absolute most boundaries of each starshipID
	# structured $boundaries{$contigID}{$starshipID}}, $hoodBegin, $hoodEnd;
	my ($boundaries) = Parse_region_boundaries($genes);
	
	# then identify captain and assign it the region boundaries
	foreach my $contigID (keys %{$genes}) {
		foreach my $regionID (keys %{$genes->{$contigID}}) {
			foreach my $geneID (keys %{$genes->{$contigID}->{$regionID}}) {
				my ($begin, $end, $strand, $tag, $annotation) = @{$genes->{$contigID}->{$regionID}->{$geneID}};
				if ($tag eq 'cap') {
					push @{$elements{$contigID}{$regionID}{$geneID}}, ${$boundaries->{$contigID}->{$regionID}}[0], ${$boundaries->{$contigID}->{$regionID}}[1], $strand, $tag;
					push @{$elements{$contigID}{$regionID}{$geneID}}, $annotation if (defined $annotation);
					$uniqueElements{$regionID} = 1;
				}
			}
		}
	}
	
	# count the number of unique elements with captains found in data
	$elementCount = scalar(keys(%uniqueElements));
	
	return(\%elements, $elementCount);
}
sub Parse_region_boundaries {
	my ($regions) = @_;
 	my %coordinates;
 	foreach my $contigID (keys %{$regions}) {
 		foreach my $regionID (keys %{$regions->{$contigID}}) {
 			my ($hoodBegin, $hoodEnd) = (100000000000000000000, 0); 
			# iterate through all coordinates, and find the upstream-most and downstream-most
 			foreach my $geneID (keys %{$regions->{$contigID}->{$regionID}}) {
 				my ($currentBegin, $currentEnd, $strand, $tag) = @{$regions->{$contigID}->{$regionID}->{$geneID}};
	 			$hoodBegin = $currentBegin if ($currentBegin < $hoodBegin);
				$hoodEnd = $currentEnd if ($currentEnd > $hoodEnd);
			}
			push @{$coordinates{$contigID}{$regionID}}, $hoodBegin, $hoodEnd;
		}
	}
	return(\%coordinates);	 
}



1;