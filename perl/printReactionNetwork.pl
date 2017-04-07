#!/usr/bin/perl -w

#################################################################
#
# This script takes a toychem reaction network output and
# generates an image of the according reaction network.
# Therein, all molecules are depicted as molecule graphs
# within their respective nodes using OpenBabel >= v2.3 for the
# generation. The network layouting is done using tools from
# the graphviz package and provided in SVG format.
# 
# author : Martin Mann (c) 2012
# 
#################################################################
#
# Call possibilities are as follows where <..> represents the 
# remaining toyChem parameters to be set by the user.
# 
#  (a) within pipeline (loss of remaining toyChem output)
#  
# $ toyChem -output=N -out=STDOUT <..> | perl printReactionNetwork.pl > NETWORK.svg
#
#  (b) via intermediate toyChem output file
# 
# $ toyChem -output=N -out=toyChem.output <..>
# $ perl printReactioNetwork.pl < toyChem.output > NETWORK.svg
#
#  (c) direct conversion to PDF via rsvg-convert
#  
# $ .. | perl printReactionNetwork.pl | rsvg-convert -f pdf > NETWORK.pdf
#
#################################################################

require File::Temp;

use strict;
use English;

use File::Temp qw/ tempdir /;


# check for graph layout program
my $graphvizFound = 0;
open (GRAPHVIZCHECK, "which neato 2>/dev/null | ");
while($graphvizFound == 0 and my $line = <GRAPHVIZCHECK>) {
	if ($line =~ m/\/neato$/) {
		$graphvizFound = 1;
	}
}
close GRAPHVIZCHECK;
if ($graphvizFound == 0) {
	die "Cannot run script: neato from graphviz package is not available";
}

# check for obabel program for molecule graph image generation
my $openbabelFound = 0;
open (OPENBABELCHECK, "which obabel 2>/dev/null | ");
while($openbabelFound == 0 and my $line = <OPENBABELCHECK>) {
if ($line =~ m/\/obabel$/) {
		$openbabelFound = 1;
	}
}
close OPENBABELCHECK;
if ($openbabelFound == 0) {
	die "Cannot run script: obabel is not available";
}

# check for rsvg-convert program for final SVG conversion
my $rsvgFound = 0;
open (RSVGCHECK, "which rsvg-convert 2>/dev/null | ");
while($rsvgFound == 0 and my $line = <RSVGCHECK>) {
	if ($line =~ m/convert$/) {
		$rsvgFound = 1;
	}
}
close (RSVGCHECK);
if ($rsvgFound == 0) {
	die "Cannot run script: rsvg-convert is not available";
}

# create temporary directory where to produce the molecule graph images
my $tempDir = tempdir( CLEANUP => 1 );
chdir $tempDir;

# data for reaction network extraction from toyChem output
my $graphStartFound = 0;
my $newGraph = "";

# read and parse toychem output from STDIN
while( my $line = <STDIN> ) {

	 # find and store graph start
	if ($graphStartFound == 0) {
		if ($line =~ m/digraph reactionNetwork/) {
			$graphStartFound = 1;
			$newGraph .= $line;
		}
		next;
	}
	
	 # handle molecule nodes
	if ($line =~ m/^\s+(M\d+)\s+\[.*label="(\S+)"\]/) {
		 # get information for graphics production
		my $molID = $1;
		my $molSMILES = $2;
		 # produce graphics
		system("obabel -:\"$molSMILES\" -O $molID.png -d 2>/dev/null;");
		 # write new node information
		$newGraph .= "  $molID [shape=oval image=\"$molID.png\" label=\"\"];\n";
	} else {
	 # handle all other graph information
		$newGraph .= $line;
	}

}

 # write new graph in DOT format
open( OUT, ">newgraph.dot");
print OUT $newGraph;
close( OUT );

 # run dot/fdp to create the picture
system("dot newgraph.dot -Tsvg -O;");

 # generate final SVG and print to stream
open( SVG, "rsvg-convert -f svg < newgraph.dot.svg |");
while (my $svg = <SVG>) {
	print STDOUT $svg;
}
close( SVG );

 # leave and destroy temporary directory

 # remove temporary data
chdir "..";


## EOF

