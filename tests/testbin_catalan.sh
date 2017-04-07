#!/bin/sh

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/catalan"
echo "==============================================";
echo

CALL=../src/bin/catalan

echo " --> checking help"
$CALL -help 2>&1

echo
echo


RUNPARAM=" -catalan=catalan/2-step.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1

RUNPARAM=" -catalan=catalan/3-step.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1


echo
echo "===============  END TEST  ===================";