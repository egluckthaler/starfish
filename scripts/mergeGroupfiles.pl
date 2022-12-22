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
-t    FILE    target MCL-formatted group file
-q    FILE    query MCL-formatted group file

Required, with defaults:
-s    STR     separator between the merged groupIDs (default: '|')\n/;

	if (not defined $message) {
		$message = qq/
This script will merge two MCL formatted groups files by their members, creating new groups
that consist of combinations of the groups present in -t and -q. Prints to STDOUT. The
groupIDs of -t will come first in the new merged groupID \n\n/; 
	} else {
		$message = "\n$message\nuse -h for more details\n" ;
	}	
	die($usage, $message);
}

main: {
	my %opts;
	getopt('t:q:s:h', \%opts);
	Opts_check(\%opts);
	
	my ($targetGroups) = Parse_group_file_by_member($opts{'t'});
	my ($queryGroups) = Parse_group_file_by_member($opts{'q'});
	
	my (%newGroups, %observedMembers);
	foreach my $targetMember (keys %{$targetGroups}) {
		if (exists $queryGroups->{$targetMember}) {
			my $newGroupID = "$targetGroups->{$targetMember}$opts{'s'}$queryGroups->{$targetMember}";
			$newGroups{$newGroupID}{$targetMember} = 1;
		} else {
			$newGroups{$targetGroups->{$targetMember}}{$targetMember} = 1;
		}
	}
	foreach my $queryMember (keys %{$queryGroups}) {
		if (not exists $targetGroups->{$queryMember}) {
			$newGroups{$queryGroups->{$queryMember}}{$queryMember} = 1;
		}
	}
	foreach my $newGroupID (sort keys %newGroups) {
		print $newGroupID;
		foreach my $memberID (sort keys %{$newGroups{$newGroupID}}) {
			print "\t$memberID";
		}
		print "\n";
	}
}

sub Parse_group_file_by_member {
	my ($clusteringOutfile) = @_;
	my $datestring = localtime();					
	my %member2group;
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cannot read $clusteringOutfile, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		my (@members) = split/\s+/, $line;
		my $group = shift @members;
		$group =~ s/:$//;
		
		foreach my $member (@members) {
			$member2group{$member} = $group;
		}
	}
	return(\%member2group);	
}


sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage(" ") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a file to -t\n") if (not defined $opts->{'t'});
	usage("\nError: the file provided to -t does not exist\n") if (! -f $opts->{'t'});
	usage("\nError: please provide a file to -q\n") if (not defined $opts->{'q'});
	usage("\nError: the file provided to -q does not exist\n") if (! -f $opts->{'q'});
	if (not defined $opts->{'s'}) {
		$opts->{'s'} = '|';
	}
}