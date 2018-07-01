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
	
	if($arg =~ /\.(fastq\w+)$/ || $arg =~ /\.(fastq\w+\.gz)$/){
		
		$format = $1;
		my $file = $arg;
		
		$arg =~ s/\.fastq\w+$/\.fastq/;
		$arg =~ s/\.fastq\w+\.gz$/\.fastq\.gz/;
		
		$ARGV[$_] = $arg;
		rename $file, $arg;
		
		push @inFiles,  $arg if $arg =~ /\.dat_input\.fastq$/ || $arg =~ /\.dat_input\.fastq\.gz$/;
		push @outFiles, $arg if $arg =~ /\.dat\.fastq$/       || $arg =~ /\.dat\.fastq\.gz$/;
	}
}

my $call = join " ", @ARGV;

system $call and exit 1;


unlink $_ or warn "Could not unlink $_: $!" foreach(@inFiles);

foreach(@outFiles){
	
	my $file = $_;
	
	s/\.fastq$//;
	s/\.fastq\.gz$//;
	
	rename $file, $_;
}

