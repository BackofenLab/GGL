#!/bin/sh

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/gameOfLife"
echo "==============================================";
echo

CALL=../src/bin/gameOfLife

echo " --> checking help"
$CALL -help 2>&1

echo
echo


RUNPARAM=" -rules=game-of-life/gof.rules -grid=game-of-life/gof.glider -iter=10"
echo " --> running glider grid"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1


echo
echo "===============  END TEST  ===================";