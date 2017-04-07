#!/usr/bin/perl

#################################################################
#
# This script takes the graph grammar GML encoding of a chemical 
# reaction in GML format and  generates an SVG depiction of the 
# according reaction transition. Therein, all 2D-embedding is 
# done via OpenBabel >= v2.3.
#
# NOTE : Ensure that the input contains only ONE RULE ENCODING!
# 
# author : Martin Mann (c) 2012
# 
#################################################################
#
# Call possibilities are as follows.
# 
#  (a) within pipeline 
#  
# $ cat RULE.gml | perl chemrule2svg.pl > RULE.svg
#
#  (b) direct call
# 
# $ perl chemrule2svg.pl < RULE.gml > RULE.svg
#
#  (c) direct conversion to pdf via 'rsvg-convert'
# 
# $ perl chemrule2svg.pl < RULE.gml | rsvg-convert -f pdf > RULE.pdf
#
#################################################################


#################################################################

use strict;
use English;
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

# arrow spacer width, the arrow length is ($arrowWidth-(2*$padding*$scaling))
my $arrowWidth = 200;


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
my $ruleID = "";
my @nodeIDLeft = ();
my @nodeIDRight = ();
my @nodeNameLeft = ();
my @nodeNameRight = ();
my @nodeLabel = ();
my @nodeContext = ();
my @edgeFrom = ();
my @edgeTo = ();
my @edgeLabel = ();
my @edgeContext = ();

my $nodeOpen = 0;
my $edgeOpen = 0;
my $leftAtomCount = 1;
my $rightAtomCount = 1;

# fill rule data
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
	}
	 # handle context change
	elsif ($item eq "left") {
		$context = 0;
	} elsif ($item eq "context") {
		$context = 1;
	} elsif ($item eq "right") {
		$context = 2;
	}
	 # handle node opening
	elsif ($item eq "node") {
		$nodeOpen = 1;
	}
	 # handle edge opening
	elsif ($item eq "edge") {
		$edgeOpen = 1;
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
			 # store data
			push( @nodeContext, $context );
			 # store node name in CML output
			if ($context != 2) {
				push( @nodeIDLeft, $data[$i] );
				push( @nodeNameLeft, "a".$leftAtomCount );
				$leftAtomCount++;
			} else {
				push( @nodeIDLeft, "na" );
				push( @nodeNameLeft, "na" );
			}
			if ($context != 0) {
				push( @nodeIDRight, $data[$i] );
				push( @nodeNameRight, "a".$rightAtomCount );
				$rightAtomCount++;
			} else {
				push( @nodeIDRight, "na" );
				push( @nodeNameRight, "na" );
			}
			$i++;
		}
	}
	 # handle edge from
	elsif ($item eq "source") {
		if ($edgeOpen == 1) {
			 # store data
			push( @edgeFrom, $data[$i] );
			push( @edgeContext, $context );
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
	elsif ($item eq "ruleID") {
         # store data
        $ruleID = $data[$i] ;
        $i++;
	}
	
	# TODO : store constraint information 
}

#print STDERR join("\t",@nodeLabel),"\n";
#print STDERR join("\t",@nodeContext),"\n";
#print STDERR join("\t",@nodeIDLeft),"\n";
#print STDERR join("\t",@nodeNameLeft),"\n";
#print STDERR join("\t",@nodeIDRight),"\n";
#print STDERR join("\t",@nodeNameRight),"\n";
#print STDERR "\n";
#print STDERR join("\t", @edgeFrom),"\n";
#print STDERR join("\t", @edgeTo ),"\n";
#print STDERR join("\t", @edgeLabel ),"\n";
#print STDERR join("\t", @edgeContext),"\n";
#print STDERR "\n";


########################################################################################  
# CREATE CML ENCODING OF LEFT AND RIGHT SIDE FOR 2D-EMBEDDING VIA OPENBABEL
########################################################################################  

# graphs to fill : init with CML start
my $left = "<?xml version=\"1.0\"?>\n<molecule xmlns=\"http://www.xml-cml.org/schema\">\n <atomArray>\n";
my $right = "<?xml version=\"1.0\"?>\n<molecule xmlns=\"http://www.xml-cml.org/schema\">\n <atomArray>\n";

my %leftID2idx = ();
my %rightID2idx = ();
my %leftName2idx = ();
my %rightName2idx = ();

# add node information
for ($i = 0; $i<=$#nodeContext; $i++) {
	 # left side node
	if ($nodeContext[$i] != 2) {
		$leftID2idx{$nodeIDLeft[$i]} = $i;
		$leftName2idx{$nodeNameLeft[$i]} = $i;
		$left .= "  <atom id=\"".$nodeNameLeft[$i]."\" elementType=\"";
		 # wildcard label handling -> replace with "C" for embedding
		if ($nodeLabel[$i] eq "*" or $nodeLabel[$i] =~ /^{/) {
			$left .= "C\" />\n";
		 # class ID handling within label
		} elsif ($nodeLabel[$i] =~ m/^(.+):\i+$/) {
			$left .= $1;
		} else {
			$left .= $nodeLabel[$i]."\" />\n";
		}
	}
	 # right side node
	if ($nodeContext[$i] != 0) {
		$rightID2idx{$nodeIDRight[$i]} = $i;
		$rightName2idx{$nodeNameRight[$i]} = $i;
		$right .= "  <atom id=\"".$nodeNameRight[$i]."\" elementType=\"";
		 # wildcard label handling -> replace with "C" for embedding
		 # group label handling -> replace with "C" for embedding
		if ($nodeLabel[$i] eq "*" || $nodeLabel[$i] =~ /^{/) {
			$right .= "C\" />\n";
		 # class ID handling within label
		} elsif ($nodeLabel[$i] =~ m/^(.+):\i+$/) {
			$left .= $1;
		} else {
			$right .= $nodeLabel[$i]."\" />\n";
		}
	}
}
# add node end, start bonds
$left .= " </atomArray>\n <bondArray>\n";
$right .= " </atomArray>\n <bondArray>\n";


# add bond information
for ($i = 0; $i<=$#edgeLabel; $i++) {
	 # left side bond
	if ($edgeContext[$i] != 2) {
		$left .= "  <bond atomRefs2=\"".$nodeNameLeft[$leftID2idx{$edgeFrom[$i]}]." ".$nodeNameLeft[$leftID2idx{$edgeTo[$i]}]."\" order=\"";
		if ($edgeLabel[$i] eq "-") {
			$left .= "1";
		} elsif ($edgeLabel[$i] eq "=") {
			$left .= "2";
		} elsif ($edgeLabel[$i] eq "#") {
			$left .= "3";
		} elsif ($edgeLabel[$i] eq ":") {
			$left .= "A";
		} else { # handling of wildcards etc.
			$left .= "1";
		}
		$left .= "\" />\n";
	}
	 # right side bond
	if ($edgeContext[$i] != 0) {
		$right .= "  <bond atomRefs2=\"".$nodeNameRight[$rightID2idx{$edgeFrom[$i]}]." ".$nodeNameRight[$rightID2idx{$edgeTo[$i]}]."\" order=\"";
		if ($edgeLabel[$i] eq "-") {
			$right .= "1";
		} elsif ($edgeLabel[$i] eq "=") {
			$right .= "2";
		} elsif ($edgeLabel[$i] eq "#") {
			$right .= "3";
		} elsif ($edgeLabel[$i] eq ":") {
			$right .= "A";
		} else { # handling of wildcards etc.
			$right .= "1";
		}
		$right .= "\" />\n";
	}
}

# add bond and CML end
$left .= " </bondArray>\n</molecule>";
$right .= " </bondArray>\n</molecule>";

# print STDERR  $left,"\n\n",$right,"\n"; 


my @leftX = ();
my @leftY = ();
my @rightX = ();
my @rightY = ();

my $minX = 99999; 
my $minY = 99999;
my $maxX = -99999;
my $maxY = -99999;


########################################################################################  
#  GET 2D-EMBEDDING VIA OPENBABEL OF LEFT AND RIGHT SIDE 
########################################################################################  

# get 2D layout via obenbabel
open ( LEFT2D, "echo '$left' | obabel -icml --gen2D -ocml -h 2>/dev/null |" );
#open ( LEFT2D, "<tmp.left.cml" );
while (my $line = <LEFT2D>) {
	if ($line =~ /<atom[^>]+id=\"(\S+)\"[^>]+\sx2=\"(\S+)\"[^>]+y2=\"(\S+)\"[^>]*>/) {
		if (exists $leftName2idx{$1}) {
			$leftX[$leftName2idx{$1}] = $2;
            $leftY[$leftName2idx{$1}] = $3;
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
while ( my ($key, $value) = each(%leftID2idx) ) {
	 # shift
     $leftX[$value] += $minX;
     $leftY[$value] += $minY;
	 # scale
     $leftX[$value] *= $scaling;
     $leftY[$value] *= $scaling;
     
     $maxX = max($leftX[$value],$maxX);
     $maxY = max($leftY[$value],$maxY);
}
$maxY += ($padding*$scaling);

my $leftWidth = $maxX;
my $leftHeight = $maxY;


$minX = 999999;
$minY = 999999;
# get 2D layout via obenbabel
open ( RIGHT2D, "echo '$right' | obabel -icml --gen2D -ocml -h 2>/dev/null |" );
#open ( RIGHT2D, "<tmp.right.cml" );
while (my $line = <RIGHT2D>) {
	if ($line =~ /<atom[^>]+id=\"(\S+)\"[^>]+\sx2=\"(\S+)\"[^>]+y2=\"(\S+)\"[^>]*>/) {
		if (exists $rightName2idx{$1}) {
            $rightX[$rightName2idx{$1}] = $2;
            $rightY[$rightName2idx{$1}] = $3;
            $minX = min($2,$minX);
            $minY = min($3,$minY);
        }
	}
}
close( RIGHT2D );

# shift all coordinate values to positive
$minX = ($minX*-1);
$minY = ($minY*-1) + $padding;
$maxX = 0;
$maxY = 0;
while ( my ($key, $value) = each(%rightID2idx) ) {
	 # shift
     $rightX[$value] += $minX;
     $rightY[$value] += $minY;
	 # scale
     $rightX[$value] *= $scaling;
     $rightY[$value] *= $scaling;
	 # shift further
     $rightX[$value] += $leftWidth + $arrowWidth;
     
     $maxX = max($rightX[$value],$maxX);
     $maxY = max($rightY[$value],$maxY);
}
$maxX += ($padding*$scaling);
$maxY += ($padding*$scaling);

my $arrowY = max($maxY,$leftHeight)/2;

$maxY = max($maxY, $leftHeight);

########################################################################################  
#  WRITE SVG OUTPUT OF THE RULE
########################################################################################  

# SVG header
$svg = "<?xml version=\"1.0\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\"\nxmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:cml=\"http://www.xml-cml.org/schema\" width=\"$maxX\" height=\"$maxY\" x=\"0\" y=\"0\" font-family=\"sans-serif\" stroke=\"rgb(0,0,0)\" stroke-width=\"1\"  stroke-linecap=\"round\">\n";

# write left side

# add bonds
for ($i = 0; $i<=$#edgeLabel; $i++) {
	 # left side bond
	if ($edgeContext[$i] != 2) {
		my $from = $leftID2idx{$edgeFrom[$i]};
		my $to = $leftID2idx{$edgeTo[$i]};
		my $x = $leftX[$from] + (($leftX[$to]-$leftX[$from])/2);
		my $y = $leftY[$from] + (($leftY[$to]-$leftY[$from])/2);
		 # write line
        $svg .= " <line x1=\"".$leftX[$from]."\" y1=\"".$leftY[$from]."\" x2=\"".$leftX[$to]."\" y2=\"".$leftY[$to]."\" stroke=\"rgb(245,245,245)\"  stroke-width=\"10\"/>\n";
         # write identifier
        $svg .= " <text x=\"".$x."\" y=\"".$y."\" ";
        if ($edgeContext[$i] == 0) {
	        $svg .= "fill=\"rgb(255,0,0)\"  stroke=\"rgb(255,0,0)\"";
        } else {
	        $svg .= "fill=\"rgb(0,0,0)\"  stroke=\"rgb(0,0,0)\"";
        }
        $svg .= " id=\"".$edgeFrom[$i]."-".$edgeTo[$i]."\"";
        $svg .= " stroke-width=\"1\" font-size=\"16\" text-anchor=\"middle\" >".$edgeLabel[$i]."</text>\n";
        
	}
}
# add nodes
while ( my ($key, $value) = each(%leftID2idx) ) {
	$svg .= " <text x=\"".$leftX[$value]."\" y=\"".$leftY[$value]."\" ";
	if ($nodeContext[$value] == 0) {
		$svg .= "fill=\"rgb(255,0,0)\"  stroke=\"rgb(255,0,0)\"";
	} else {
		$svg .= "fill=\"rgb(0,0,0)\"  stroke=\"rgb(0,0,0)\"";
	}
	$svg .= " stroke-width=\"1\" font-size=\"16\" text-anchor=\"middle\" >";
	$svg .= "".$nodeLabel[$value];
	$svg .= "<tspan font-size=\"10\" stroke-width=\"0.7\" >".$nodeIDLeft[$value]."</tspan>\n";
	$svg .=" </text>\n";
}

# write arrow

$svg .= " <polygon fill=\"rgb(0,0,0)\" points=\"";
$svg .= ($leftWidth+($padding*$scaling)).",".($arrowY-2)." ";
$svg .= ($leftWidth+$arrowWidth-($padding*$scaling)-20).",".($arrowY-2)." ";
$svg .= ($leftWidth+$arrowWidth-($padding*$scaling)-20).",".($arrowY-7)." ";
$svg .= ($leftWidth+$arrowWidth-($padding*$scaling)).",".($arrowY)." ";
$svg .= ($leftWidth+$arrowWidth-($padding*$scaling)-20).",".($arrowY+7)." ";
$svg .= ($leftWidth+$arrowWidth-($padding*$scaling)-20).",".($arrowY+2)." ";
$svg .= ($leftWidth+($padding*$scaling)).",".($arrowY+2)." ";
$svg .= ($leftWidth+($padding*$scaling)).",".($arrowY-2)." ";
$svg .= "\"/>";

# write right side

# add bonds
for ($i = 0; $i<=$#edgeLabel; $i++) {
	 # left side bond
	if ($edgeContext[$i] != 0) {
		my $from = $rightID2idx{$edgeFrom[$i]};
		my $to = $rightID2idx{$edgeTo[$i]};
		my $x = $rightX[$from] + (($rightX[$to]-$rightX[$from])/2);
		my $y = $rightY[$from] + (($rightY[$to]-$rightY[$from])/2);
		 # write line
        $svg .= " <line x1=\"".$rightX[$from]."\" y1=\"".$rightY[$from]."\" x2=\"".$rightX[$to]."\" y2=\"".$rightY[$to]."\" stroke=\"rgb(245,245,245)\"  stroke-width=\"10\"/>\n";
         # write identifier
        $svg .= " <text x=\"".$x."\" y=\"".$y."\" ";
        if ($edgeContext[$i] == 2) {
	        $svg .= "fill=\"rgb(50,205,50)\"  stroke=\"rgb(50,205,50)\"";
        } else {
	        $svg .= "fill=\"rgb(0,0,0)\"  stroke=\"rgb(0,0,0)\"";
        }
        $svg .= " id=\"".$edgeFrom[$i]."-".$edgeTo[$i]."\"";
        $svg .= " stroke-width=\"1\" font-size=\"16\" text-anchor=\"middle\" >".$edgeLabel[$i]."</text>\n";
	}
}
# add nodes
while ( my ($key, $value) = each(%rightID2idx) ) {
	$svg .= " <text x=\"".$rightX[$value]."\" y=\"".$rightY[$value]."\" ";
	if ($nodeContext[$value] == 2) {
		$svg .= "fill=\"rgb(50,205,50)\"  stroke=\"rgb(50,205,50)\"";
	} else {
		$svg .= "fill=\"rgb(0,0,0)\"  stroke=\"rgb(0,0,0)\"";
	}
	$svg .= " stroke-width=\"1\" font-size=\"16\" text-anchor=\"middle\" >";
	$svg .= "".$nodeLabel[$value];
	$svg .= "<tspan font-size=\"10\" stroke-width=\"0.7\" >".$nodeIDRight[$value]."</tspan>\n";
	$svg .=" </text>\n";
}

# SVG footer
$svg .= "<title>$ruleID</title>\n</svg>\n";

# write SVG output to stream
print STDOUT $svg;

## EOF


