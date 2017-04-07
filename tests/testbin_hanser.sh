#!/bin/sh

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/hanser"
echo "==============================================";
echo

CALL=../src/bin/hanser

echo " --> checking help"
$CALL -help 2>&1

echo
echo


RUNPARAM=" -graph=hanser/single-loop.input -v"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1

RUNPARAM=" -graph=hanser/two-fused-rings.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1

RUNPARAM=" -graph=hanser/three-fused-rings.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1


echo
echo "===============  END TEST  ===================";