#!/usr/bin/perl

#################################################################
#
# This script takes a molecule componentn in GML encoding 
# and generates an according SVG depiction. 
# Therein, all 2D-embedding is done via OpenBabel >= v2.3.
#
# NOTE : Ensure that the input contains only ONE COMPONENT ENCODING!
# 
# author : Martin Mann (c) 2012
# 
#################################################################
#
# Call possibilities are as follows.
# 
#  (a) within pipeline 
#  
# $ cat GROUP.gml | perl molcomp2svg.pl > GROUPS.svg
#
#  (b) direct call
# 
# $ perl molcomp2svg.pl < GROUPS.gml > GROUPS.svg
#
#  (c) direct conversion to pdf via 'rsvg-convert'
# 
# $ molcomp2svg.pl < GROUPS.gml | perl rsvg-convert -f pdf > GROUPS.pdf
#
#################################################################


#################################################################

use strict;
use List::Util qw[min max];


# check for obabel program for 2D-embedding
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


#################################################################
#  SVG CREATION PARAMETERS
#################################################################

# scaling factor for the 2D-embedding coordinates of OpenBabel in CML-format
my $scaling = 50;

# relative padding according to the 2D-embedding coordinates of OpenBabel in CML-format
my $padding = 1;


#################################################################
#  INPUT PARSING
#################################################################


# container to temporary store input data
my @input = ();

# read and parse GML input from STDIN
while( my $line = <STDIN> ) {

	# ensure spaces before and after brackets
	$line =~ s/\[/ \[ /g;
	$line =~ s/\]/ \] /g;

	if ($line =~ /^\s*(\S.*\S)\s*$/) {
		push ( @input, split( /[ \t\n]+/, $1 ) );
	}

}


########################################################################################  
#  PARSE GML RULE ENCODING AND STORE INFORMATION
########################################################################################  

 # container that will hold all entries
my @data = ();
 # test and copy all GML entries
my $i=0;
my $n=$i;
while($i <= $#input) {
	 # handle incomplete string labels
	if ($input[$i] =~ /^["'].*[^"']$/) {
		 # add next element to this string label
		$input[$i] .= " $input[$n]";
		 # remove added label part
		$input[$n] = "";
		 # increase counter of label part
		$n++;
	} # prune " of string labels
	elsif ($input[$i] =~ /^["'](.*)["']$/) {
		$input[$i] = $1;
	} else {
		 # copy non-empty entries
		if ($input[$i] =~ /\S+/) {
			push( @data, $input[$i] );
		}
		 # go to next item within input
		$i++;
		$n = $i+1;
	}
}
# clear input data container (not needed anymore
@input = ();

# context encoding : left=0, context=1, right=2
my $context = 0;

# print STDERR join ("\n",@data),"\n\n";

# parsed data
my $description = "";
my @nodeID = ();
my @nodeName = ();
my @nodeLabel = ();
my @edgeFrom = ();
my @edgeTo = ();
my @edgeLabel = ();
my @compIDs = ();
my @nodeIsCompID = ();

my $nodeOpen = 0;
my $edgeOpen = 0;
my $compIDOpen = 0;
my $atomCount = 1;

# fill molComp data
$i = 0;
while ($i < $#data) {
	 # access to data element
	my $item = $data[$i];
	$i++;
	
	 # ignore opening brackets
	if ($item =~ /^[\[]$/) {
		next;
	}
	 # handle closing brackets
	elsif ($item eq "]") {
		if ($nodeOpen == 1) {
			$nodeOpen = 0;
		}
		if ($edgeOpen == 1) {
			$edgeOpen = 0;
		}
		if ($compIDOpen == 1) {
			$compIDOpen = 0;
		}
	}
	 # handle node opening
	elsif ($item eq "node") {
		$nodeOpen = 1;
	}
	 # handle edge opening
	elsif ($item eq "edge") {
		$edgeOpen = 1;
	}
	 # handle compIDs opening
	elsif ($item eq "compIDs") {
		$compIDOpen = 1;
	}
	 # handle label
	elsif ($item eq "label") {
		if ($nodeOpen == 1) {
			 # store data
			push( @nodeLabel, $data[$i] );
			$i++;
		} elsif ($edgeOpen == 1) {
			 # store data
			push( @edgeLabel, $data[$i] );
			$i++;
		}
	}
	 # handle node id
	elsif ($item eq "id") {
		if ($nodeOpen == 1) {
			 # store node name in CML output
			push( @nodeID, $data[$i] );
			push( @nodeIsCompID, 0 );
			push( @nodeName, "a".$atomCount );
			$atomCount++;
			$i++;
		}
		if ($compIDOpen == 1) {
			push( @compIDs, $data[$i] );
			$i++;
		}
	}
	 # handle edge from
	elsif ($item eq "source") {
		if ($edgeOpen == 1) {
			 # store data
			push( @edgeFrom, $data[$i] );
			$i++;
		}
	}
	 # handle edge to
	elsif ($item eq "target") {
		if ($edgeOpen == 1) {
			 # store data
			push( @edgeTo, $data[$i] );
			$i++;
		}
	}
	 # ruleID
	elsif ($item eq "description") {
         # store data
        $description = $data[$i] ;
        $i++;
	}
	
	# TODO : store constraint information 
}

#print STDERR join("\t",@nodeLabel),"\n";
#print STDERR join("\t",@nodeID),"\n";
#print STDERR join("\t",@nodeName),"\n";
#print STDERR "\n";
#print STDERR join("\t", @edgeFrom),"\n";
#print STDERR join("\t", @edgeTo ),"\n";
#print STDERR join("\t", @edgeLabel ),"\n";
#print STDERR "\n";
#print STDERR join("\t", @compIDs ),"\n";
#print STDERR "\n";


########################################################################################  
# CREATE CML ENCODING FOR 2D-EMBEDDING VIA OPENBABEL
########################################################################################  

# graphs to fill : init with CML start
my $cml = "<?xml version=\"1.0\"?>\n<molecule xmlns=\"http://www.xml-cml.org/schema\">\n <atomArray>\n";

my %cmlID2idx = ();
my %cmlName2idx = ();

# add node information
for ($i = 0; $i<=$#nodeID; $i++) {
	$cmlID2idx{$nodeID[$i]} = $i;
	$cmlName2idx{$nodeName[$i]} = $i;
	$cml .= "  <atom id=\"".$nodeName[$i]."\" elementType=\"";
	 # wildcard label handling -> replace with "C" for embedding
	if ($nodeLabel[$i] eq "*" or $nodeLabel[$i] =~ /^{/) {
		$cml .= "C\" />\n";
	 # class ID handling within label
	} elsif ($nodeLabel[$i] =~ m/^(.+):\i+$/) {
		$cml .= $1;
	} else {
		$cml .= $nodeLabel[$i]."\" />\n";
	}
}
# add node end, start bonds
$cml .= " </atomArray>\n <bondArray>\n";

# add bond information
for ($i = 0; $i<=$#edgeLabel; $i++) {
	$cml .= "  <bond atomRefs2=\"".$nodeName[$cmlID2idx{$edgeFrom[$i]}]." ".$nodeName[$cmlID2idx{$edgeTo[$i]}]."\" order=\"";
	if ($edgeLabel[$i] eq "-") {
		$cml .= "1";
	} elsif ($edgeLabel[$i] eq "=") {
		$cml .= "2";
	} elsif ($edgeLabel[$i] eq "#") {
		$cml .= "3";
	} elsif ($edgeLabel[$i] eq ":") {
		$cml .= "A";
	} else { # handling of wildcards etc.
		$cml .= "1";
	}
	$cml .= "\" />\n";
}

# add bond and CML end
$cml .= " </bondArray>\n</molecule>";

# print STDERR  $cml,"\n"; 


my @cmlX = ();
my @cmlY = ();

my $minX = 99999; 
my $minY = 99999;
my $maxX = -99999;
my $maxY = -99999;


########################################################################################  
#  GET 2D-EMBEDDING VIA OPENBABEL OF LEFT AND RIGHT SIDE 
########################################################################################  

# get 2D layout via obenbabel
open ( LEFT2D, "echo '$cml' | obabel -icml --gen2D -ocml -h 2>/dev/null |" );
#open ( LEFT2D, "<tmp.left.cml" );
while (my $line = <LEFT2D>) {
	if ($line =~ /<atom[^>]+id=\"(\S+)\"[^>]+\sx2=\"(\S+)\"[^>]+y2=\"(\S+)\"[^>]*>/) {
		if (exists $cmlName2idx{$1}) {
			$cmlX[$cmlName2idx{$1}] = $2;
            $cmlY[$cmlName2idx{$1}] = $3;
            $minX = min($2,$minX);
            $minY = min($3,$minY);
		}
	}
}
close( LEFT2D );

# generate SVG output with coordinates
my $svg = "";

# shift all coordinate values to positive
$minX = ($minX*-1) + $padding;
$minY = ($minY*-1) + $padding;
$maxX = 0;
$maxY = 0;
while ( my ($key, $value) = each(%cmlID2idx) ) {
	 # shift
     $cmlX[$value] += $minX;
     $cmlY[$value] += $minY;
	 # scale
     $cmlX[$value] *= $scaling;
     $cmlY[$value] *= $scaling;
     
     $maxX = max($cmlX[$value],$maxX);
     $maxY = max($cmlY[$value],$maxY);
}
$maxX += ($padding*$scaling);
$maxY += ($padding*$scaling);


########################################################################################  
#  GET COMPID INFORMATION
########################################################################################  

# mark all compID nodes
foreach my $id (@compIDs) {
	$nodeIsCompID[ $cmlID2idx{$id} ] = 1;
}

########################################################################################  
#  WRITE SVG OUTPUT OF THE MOLCOMP
########################################################################################  

# SVG header
$svg = "<?xml version=\"1.0\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\"\nxmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:cml=\"http://www.xml-cml.org/schema\" width=\"$maxX\" height=\"$maxY\" x=\"0\" y=\"0\" font-family=\"sans-serif\" stroke=\"rgb(0,0,0)\" stroke-width=\"1\"  stroke-linecap=\"round\">\n";

# write left side

# add bonds
for ($i = 0; $i<=$#edgeLabel; $i++) {
		my $from = $cmlID2idx{$edgeFrom[$i]};
		my $to = $cmlID2idx{$edgeTo[$i]};
		my $x = $cmlX[$from] + (($cmlX[$to]-$cmlX[$from])/2);
		my $y = $cmlY[$from] + (($cmlY[$to]-$cmlY[$from])/2);
		 # write line
        $svg .= " <line x1=\"".$cmlX[$from]."\" y1=\"".$cmlY[$from]."\" x2=\"".$cmlX[$to]."\" y2=\"".$cmlY[$to]."\" stroke=\"rgb(245,245,245)\"  stroke-width=\"10\"/>\n";
         # write identifier
        $svg .= " <text x=\"".$x."\" y=\"".$y."\" ";
        $svg .= "fill=\"rgb(0,0,0)\"  stroke=\"rgb(0,0,0)\"";
        $svg .= " id=\"".$edgeFrom[$i]."-".$edgeTo[$i]."\"";
        $svg .= " stroke-width=\"1\" font-size=\"16\" text-anchor=\"middle\" >".$edgeLabel[$i]."</text>\n";
        
}
# add nodes
while ( my ($key, $value) = each(%cmlID2idx) ) {
	$svg .= " <text x=\"".$cmlX[$value]."\" y=\"".$cmlY[$value]."\" ";
	if ($nodeIsCompID[$value] != 0) {
		$svg .= "fill=\"rgb(255,0,0)\"  stroke=\"rgb(255,0,0)\"";
	} else {
		$svg .= "fill=\"rgb(0,0,0)\"  stroke=\"rgb(0,0,0)\"";
	}
	$svg .= " stroke-width=\"1\" font-size=\"16\" text-anchor=\"middle\" >";
	$svg .= "".$nodeLabel[$value];
	$svg .= "<tspan font-size=\"10\" stroke-width=\"0.7\" >".$nodeID[$value]."</tspan>\n";
	$svg .=" </text>\n";
}

# SVG footer
$svg .= "<title>$description</title>\n</svg>\n";

# write SVG output to stream
print STDOUT $svg;

## EOF


