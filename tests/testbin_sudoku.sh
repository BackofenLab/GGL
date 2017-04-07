#!/bin/sh

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/sudoku"
echo "==============================================";
echo

CALL=../src/bin/sudoku

echo " --> checking help"
$CALL -help 2>&1

echo
echo


RUNPARAM=" -sudoku=sudoku/very-easy.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1

RUNPARAM=" -sudoku=sudoku/easy.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1

RUNPARAM=" -sudoku=sudoku/medium.input"
echo " --> running solver"
echo " --> parameter = $RUNPARAM"
$CALL $RUNPARAM 2>&1


echo
echo "===============  END TEST  ===================";