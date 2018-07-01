#!/usr/bin/env perl

# Flexbar wrapper for Galaxy tool definition, version 3.4.2
# Author: Johannes Roehr

use warnings;
use strict;

my $format;
my @inFiles;
my @outFiles;

foreach(0..$#ARGV){
	my $arg = $ARGV[$_];
	
	if($arg =~ /\.(fastq\w+)$/ || $arg =~ /\.(fastq\w+\.gz)$/ || $arg =~ /\.(fastq\w+\.bz2)$/){
		
		if(defined $format && $format ne $1){
			print STDERR "Paired read files should have the same format.\n";
			exit 1;
		}
		$format = $1;
		
		my $file = $arg;
		
		$arg =~ s/\.fastq\w+$/\.fastq/;
		$arg =~ s/\.fastq\w+\.gz$/\.fastq\.gz/;
		$arg =~ s/\.fastq\w+\.bz2$/\.fastq\.bz2/;
		
		push @inFiles,  $arg if $arg =~ /\.dat_input\.fastq$/ || $arg =~ /\.dat_input\.fastq\.gz$/ || $arg =~ /\.dat_input\.fastq\.bz2$/;
		push @outFiles, $arg if $arg =~ /\.dat\.fastq$/       || $arg =~ /\.dat\.fastq\.gz$/       || $arg =~ /\.dat\.fastq\.bz2$/;
		
		$ARGV[$_] = $arg;
		rename $file, $arg;
	}
}

my $call = join " ", @ARGV;

system $call and exit 1;


unlink $_ or warn "Could not unlink $_: $!" foreach(@inFiles);

foreach(@outFiles){
	
	my $file = $_;
	
	s/\.fastq$//;
	s/\.fastq\.gz$//;
	s/\.fastq\.bz2$//;
	
	rename $file, $_;
}

