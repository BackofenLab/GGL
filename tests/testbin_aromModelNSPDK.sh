#!/bin/sh

CALL_SUFFIX=" 2>&1";

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/aromModelNSPDK"
echo "==============================================";
echo

BINARY=../src/bin/aromModelNSPDK

echo " --> checking help"
$BINARY -help 2>&1

echo
echo

# apply OpenBabel model
PRECALL="cat toyChem/purine.smi | grep -E "^#" -v | "
RUNPARAM="-mode=apply -model=O"
POSTCALL=" "
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"

## aromatic molecule node weights
#PRECALL=" echo 'c1c2c(nc[nH]2)ncn1' | "
#RUNPARAM=" -model=M -mode=apply -outNodeWeight"
#POSTCALL=" | sed 's/]/]\n/g' "
#FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
#echo
#echo " --> precall   = $PRECALL"
#echo " --> parameter = $RUNPARAM"
#echo " --> postcall  = $POSTCALL"
#echo " --> full call = $FULLCALL"
#echo
#bash -c "$FULLCALL"
#
## non-aromatic molecule node weights
#PRECALL=" echo 'C1CCCC#CCC1' | "
#RUNPARAM=" -model=M -mode=apply -outNodeWeight"
#POSTCALL=" | sed 's/]/]\n/g' "
#FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
#echo
#echo " --> precall   = $PRECALL"
#echo " --> parameter = $RUNPARAM"
#echo " --> postcall  = $POSTCALL"
#echo " --> full call = $FULLCALL"
#echo
#bash -c "$FULLCALL"

# problematic molecule node weights
PRECALL=" echo 'c1nc2c(=O)[nH]c(nc2n1[CH]3[CH]([CH]([CH](O3)CO)O)O)N' | "
RUNPARAM=" -model=M -mode=apply -outNodeWeight"
POSTCALL=" | sed 's/]/]\n/g' "
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"


echo
echo "===============  END TEST  ===================";