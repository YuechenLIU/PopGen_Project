#! /usr/bin/perl

## A script to calculate the amount of
## missing data in an aligned sequence
## matrix.
## S. Hedtke  Oct 2011

## Usage : perl calculatemissing.pl filename
## Your file can be in nexus or phylip format.

my $size=0;
my $qs=0;

while (<ARGV>) {
	my $t=$_;
	chomp $t;
	if ($t eq /\s/ || $t=~/^\#/ || $t=~/\;/) {next;}
	my @poop=split(/\s/,$t);
	my $thingy=scalar@poop-1;
	if ($thingy<1) {next;}
	my $seq=$poop[$thingy];
	my @seqs=split('',$seq);
	for (my $i=0; $i<scalar@seqs; $i++) {
	$size++;
	if ($seqs[$i]=~/(\?|\-|N)/i) { $qs++; }	
	}
}
my $percent=100*$qs/$size;
print "\n Size of matrix: $size\n Number of empty cells, gaps, or Ns: $qs\nPercent missing = $percent\n";

exit;

#  Copyright (C) 2011  Shannon M Hedtke
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   See <http://www.gnu.org/licenses/>.