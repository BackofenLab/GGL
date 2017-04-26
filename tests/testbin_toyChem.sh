#!/bin/sh

echo "==============  BEGIN TEST  ==================";
echo "==============================================";
echo "               bin/toyChem"
echo "==============================================";
echo

TOYCHEM=../src/bin/toyChem

echo " --> checking help"
$TOYCHEM -help 2>&1

echo
echo


RUNPARAM=" -smiles=dielsalder/da1.smi:dielsalder/da2.smi -rules=dielsalder/da_rule.gml -iter=2 -outMode=R -v"
echo " --> running Diels-Alder reaction"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1

RUNPARAM=" -mols=dielsalder/da1.gml:dielsalder/da2.gml -rules=dielsalder/da_rule.gml -iter=2 -outMode=R -rate=A -v"
echo " --> running Diels-Alder reaction"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1

RUNPARAM=" -mols=dielsalder/da1.gml:dielsalder/da2.gml -rules=dielsalder/da_rule.gml -iter=2 -outMode=N -rate=A -v"
echo " --> running Diels-Alder reaction"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1

RUNPARAM=" -mols=dielsalder/da1.gml:dielsalder/da2.gml -rules=dielsalder/da_rule.gml -iter=2 -outMode=S -v"
echo " --> running Diels-Alder reaction"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1

RUNPARAM=" -smiles=dielsalder/da3.smi:dielsalder/da2.smi -rules=dielsalder/da_rule.gml -iter=2 -outMode=S -v"
echo " --> running Diels-Alder reaction with atom class label without atom-class-ignore"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1

RUNPARAM=" -smiles=dielsalder/da3.smi:dielsalder/da2.smi -rules=dielsalder/da_rule.gml -iter=2 -outMode=S -ignoreAtomClass -v"
echo " --> running Diels-Alder reaction with atom class label with atom-class-ignore"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1


RUNPARAM=" -smiles=toyChem/NADH.smi:toyChem/lactat.smi -groups=toyChem/groups.gml -rules=toyChem/lactat-dehydrogenase.gml -iter=1"
echo " --> running Lactat-dehydrogenase reaction"
echo " --> parameter = $RUNPARAM"
$TOYCHEM $RUNPARAM 2>&1


echo
echo "===============  END TEST  ===================";