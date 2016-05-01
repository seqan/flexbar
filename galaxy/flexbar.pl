#!/usr/bin/env perl

# Flexbar wrapper for Galaxy tool definition, version 2.5
# Author: Johannes Roehr

use warnings;
use strict;

my ($outFile, $id, $folder, $format) = @ARGV[($#ARGV - 3) .. $#ARGV];

my $call = join " ", @ARGV[0..($#ARGV - 4)];

system $call .' --target FlexbarTargetFile > '. $outFile and exit 1;


foreach(<FlexbarTargetFile*>){
	
	my $fileType;
	
	$fileType = $1         if /\.(\w+)$/;
	$fileType = $format    if /\.\w*fast\w$/;
	$fileType = 'fasta'    if /\.fasta$/;
	$fileType = 'csfasta'  if /\.csfasta$/;
	$fileType = 'tabular'  if /\.lengthdist$/;
	
	my $file = $_;
	
	s/_/-/g;
	
	my $name = "primary_". $id ."_". $_ ."_visible_". $fileType;
	
	rename $file, $name;
	rename $name, $folder;
}

