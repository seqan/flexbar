#!/usr/bin/env perl

# Flexbar wrapper for Galaxy tool definition, version 3.5.0
# Author: Johannes Roehr

use warnings;
use strict;

my $format;
my @inFiles;
my @outFiles;

my $compression = "";

foreach(0..$#ARGV){
	my $arg = $ARGV[$_];
	
	if($arg =~ /\.(fastq\w+)$/ || $arg =~ /\.(fastq\w+\.gz)$/ || $arg =~ /\.(fastq\w+\.bz2)$/){
		
		if(defined $format && $format ne $1){
			warn "Read files should have the same format.\n";
			exit 1;
		}
		$format = $1;
		
		my $file = $arg;
		
		$arg =~ s/\.fastq\w+$/\.fastq/;
		$arg =~ s/\.fastq\w+\.gz$/\.fastq\.gz/;
		$arg =~ s/\.fastq\w+\.bz2$/\.fastq\.bz2/;
		
		$compression = "GZ"  if $arg =~ /\.fastq\.gz$/;
		$compression = "BZ2" if $arg =~ /\.fastq\.bz2$/;
		
		$ARGV[$_] = $arg;
		
		if($arg =~ /\.dat_input\w\.fastq$/ || $arg =~ /\.dat_input\w\.fastq\.gz$/ || $arg =~ /\.dat_input\w\.fastq\.bz2$/){
			push @inFiles, $arg;
			rename $file, $arg;
		}
		
		push @outFiles, $arg if $arg =~ /\.dat\.fastq$/ || $arg =~ /\.dat\.fastq\.gz$/ || $arg =~ /\.dat\.fastq\.bz2$/;
	}
}

my $barcoded = 0;

$barcoded = 1 if $ARGV[$#ARGV] =~ /barcoded$/;

my $call = join " ", @ARGV[0..($#ARGV - $barcoded)];

# $call = $call ." --zip-output ". $compression if $barcoded && $compression ne "";

system $call and exit 1;


unlink $_ or warn "Could not unlink $_: $!" foreach(@inFiles);

if($barcoded){
	$format =~ s/\.gz//;
	$format =~ s/\.bz2//;
	
	foreach(<$ARGV[$#ARGV]/flexbarOut*.fastq*>){
		
		my $file = $_;
		
		s/fastq$/$format/;
		# s/fastq\.gz$/$format/;
		# s/fastq\.bz2$/$format/;
		
		rename $file, $_;
	}
}
else{
	foreach(@outFiles){
		
		my $file = $_;
		
		s/\.fastq$//;
		s/\.fastq\.gz$//;
		s/\.fastq\.bz2$//;
		
		rename $file, $_;
	}
}

