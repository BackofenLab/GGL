#ifndef DATARULE_1_ICC_
#define DATARULE_1_ICC_

#include <string>

const std::string ruleDielsAlderGML =
  "rule [\n\
	 ruleID \"Diels-Alder reaction\"\n\
	 context [\n\
	   node [ id 0 label \"C\" ]\n\
	   node [ id 1 label \"C\" ]\n\
	   node [ id 2 label \"C\" ]\n\
	   node [ id 3 label \"C\" ]\n\
	   node [ id 4 label \"C\" ]\n\
	   node [ id 5 label \"C\" ]\n\
	 ]\n\
	 left [\n\
	   edge [ source 0 target 1 label \"=\" ]\n\
	   edge [ source 1 target 2 label \"-\" ]\n\
	   edge [ source 2 target 3 label \"=\" ]\n\
	   edge [ source 4 target 5 label \"=\" ]\n\
       constrainNoEdge [ source 3 target 5 ]\n\
       constrainNoEdge [ source 0 target 4 ]\n\
	 ]\n\
	 right [\n\
	   edge [ source 0 target 1 label \"-\" ]\n\
	   edge [ source 0 target 5 label \"-\" ]\n\
	   edge [ source 1 target 2 label \"=\" ]\n\
	   edge [ source 2 target 3 label \"-\" ]\n\
	   edge [ source 3 target 4 label \"-\" ]\n\
	   edge [ source 4 target 5 label \"-\" ]\n\
	 ]\n\
	]\n\
	";

#endif /*DATARULE_1_ICC_*/
