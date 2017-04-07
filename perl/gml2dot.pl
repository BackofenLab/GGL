#!/usr/bin/perl


#################################################################
#
# This script takes a graph in GML format and converts it to DOT
# format. The latter can be visualized via the graphviz package.
# 
# author : Martin Mann (c) 2012
# 
#################################################################
#
# Call possibilities are as follows.
# 
#  (a) within pipeline 
#  
# $ cat GRAPH.gml | perl gml2dot.pl > GRAPH.dot
#
#  (b) via file input
# 
# $ perl gml2dot.pl < GRAPH.gml > GRAPH.dot
#
#  (c) direct conversion to PNG via dot from graphviz
#  
# $ .. | perl gml2dot.pl | dot -Tpng -oGRAPH.png
#
#  (d) direct conversion to SVG via dot from graphviz
#  
# $ .. | perl gml2dot.pl | dot -Tsvg -oGRAPH.svg
#
#################################################################

use strict;


# read
my $gml = "";
while (my $line = <STDIN>) {
	$gml .= $line;
}

# DOT header
print STDOUT "graph graphFromGML {\n";
# nodes
while ($gml =~ /node\s*\[(.+)\]/g) {
	my $node = $1;
	# get node ID
	if ($node =~ /id\s+(\d+)\s+/) {
		my $nodeID = $1;
		# get node label
		if ($node =~ /label\s"(.*)"/) {
			my $nodeLabel = $1;
			
			# print node in DOT
			print STDOUT "  N$nodeID [shape=oval label=\"$nodeLabel\"];\n"
		}
	}
}
# edges
while ($gml =~ /edge\s*\[(.+)\]/g) {
	my $edge = $1;
	# get edge source
	if ($edge =~ /source\s+(\d+)\s+/) {
		my $source = $1;
		if ($edge =~ /target\s+(\d+)\s+/) {
			my $target = $1;
			# get edge label
			if ($edge =~ /label\s"(.*)"/) {
				my $label = $1;
				
				# print node in DOT
				print STDOUT "  N$source -- N$target [label=\"$label\"];\n"
			}
		}
	}
}
# DOT footer
print STDOUT "} // visualization: eg. ' | dot -Tpng -O GRAPHOUTPUTFILE'"

# EOF