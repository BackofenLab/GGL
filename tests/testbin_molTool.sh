#!/bin/sh

CALL_SUFFIX=" 2>&1";

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/molTool"
echo "==============================================";
echo

BINARY=../src/bin/molTool

echo " --> checking help"
$BINARY -help 2>&1

echo
echo


PRECALL="cat toyChem/purine.smi | grep -E "^#" -v | "
RUNPARAM="-aromaticity=O"
POSTCALL=" | env LC_ALL=C sort -d -s"
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"

PRECALL="cat toyChem/purine.smi | grep -E "^#" -v | "
RUNPARAM="-aromaticity=M"
POSTCALL=" | env LC_ALL=C sort -d -s"
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"

PRECALL="head -n 2 toyChem/purine.smi | tail -n 1 | "
RUNPARAM="-outMode=G"
POSTCALL=""
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"

PRECALL="head -n 3 toyChem/purine.smi | tail -n 1 | "
RUNPARAM="-outMode=G"
POSTCALL=""
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"

PRECALL="head -n 10 toyChem/purine.smi | grep -E "^#" -v | $BINARY -outMode=G | "
RUNPARAM="-inMode=G -outMode=S"
POSTCALL=" | env LC_ALL=C sort -d -s"
FULLCALL="$PRECALL $BINARY $RUNPARAM $POSTCALL 2>&1"
echo
echo " --> precall   = $PRECALL"
echo " --> parameter = $RUNPARAM"
echo " --> postcall  = $POSTCALL"
echo " --> full call = $FULLCALL"
echo
bash -c "$FULLCALL"

PRECALL=""
RUNPARAM="-in=toyChem/NADH.smi -groups=toyChem/groups.gml -noGroupInsertion"
POSTCALL=" | env LC_ALL=C sort -d -s"
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