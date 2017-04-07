
#include <iostream>
#include <sstream>

#include "ggl/chem/MoleculeComponent.hh"
#include "ggl/chem/MoleculeComponent_GMLparser.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILESwriter.hh"
#include "ggl/chem/MoleculeUtil.hh"

#include "sgm/Graph_boost.hh"

#include "utilPrintGraph_Interface.icc"

using namespace ggl;
using namespace chem;

const std::string gml[] = {
"molcomp ["
" description \"{CONH2}\""
" compIDs [ id 0 ]"
" node [ id 0 label \"C\" ]"
" node [ id 1 label \"O\" ]"
" node [ id 2 label \"N\" ]"
" node [ id 3 label \"H\" ]"
" node [ id 4 label \"H\" ]"
" edge [ source 0 target 1 label \"=\" ]"
" edge [ source 0 target 2 label \"-\" ]"
" edge [ source 2 target 3 label \"-\" ]"
" edge [ source 2 target 4 label \"-\" ]"
"]"
,
"molcomp ["
" description \"{Ribo-ADP}\""
" compIDs [ id 0 ]"
" node [ id 0 label \"C\" ]"
" node [ id 1 label \"O\" ]"
" node [ id 2 label \"C\" ]"
" node [ id 3 label \"C\" ]"
" node [ id 4 label \"O\" ]"
" node [ id 5 label \"P\" ]"
" node [ id 6 label \"O\" ]"
" node [ id 7 label \"O\" ]"
" node [ id 8 label \"O\" ]"
" node [ id 9 label \"P\" ]"
" node [ id 10 label \"O\" ]"
" node [ id 11 label \"O\" ]"
" node [ id 12 label \"O\" ]"
" node [ id 13 label \"C\" ]"
" node [ id 14 label \"C\" ]"
" node [ id 15 label \"O\" ]"
" node [ id 16 label \"C\" ]"
" node [ id 17 label \"C\" ]"
" node [ id 18 label \"O\" ]"
" node [ id 19 label \"C\" ]"
" node [ id 20 label \"O\" ]"
" node [ id 21 label \"n\" ]"
" node [ id 22 label \"c\" ]"
" node [ id 23 label \"n\" ]"
" node [ id 24 label \"c\" ]"
" node [ id 25 label \"c\" ]"
" node [ id 26 label \"N\" ]"
" node [ id 27 label \"n\" ]"
" node [ id 28 label \"c\" ]"
" node [ id 29 label \"n\" ]"
" node [ id 30 label \"c\" ]"
" node [ id 31 label \"C\" ]"
" node [ id 32 label \"O\" ]"
" node [ id 33 label \"C\" ]"
" node [ id 34 label \"O\" ]"
" node [ id 35 label \"H\" ]"
" node [ id 36 label \"H\" ]"
" node [ id 37 label \"H\" ]"
" node [ id 38 label \"H\" ]"
" node [ id 39 label \"H\" ]"
" node [ id 40 label \"H\" ]"
" node [ id 41 label \"H\" ]"
" node [ id 42 label \"H\" ]"
" node [ id 43 label \"H\" ]"
" node [ id 44 label \"H\" ]"
" node [ id 45 label \"H\" ]"
" node [ id 46 label \"H\" ]"
" node [ id 47 label \"H\" ]"
" node [ id 48 label \"H\" ]"
" node [ id 49 label \"H\" ]"
" node [ id 50 label \"H\" ]"
" node [ id 51 label \"H\" ]"
" node [ id 52 label \"H\" ]"
" node [ id 53 label \"H\" ]"
" node [ id 54 label \"H\" ]"
" node [ id 55 label \"H\" ]"
" edge [ source 1 target 0 label \"-\" ]"
" edge [ source 2 target 1 label \"-\" ]"
" edge [ source 3 target 2 label \"-\" ]"
" edge [ source 4 target 3 label \"-\" ]"
" edge [ source 5 target 4 label \"-\" ]"
" edge [ source 6 target 5 label \"-\" ]"
" edge [ source 7 target 5 label \"=\" ]"
" edge [ source 8 target 5 label \"-\" ]"
" edge [ source 9 target 8 label \"-\" ]"
" edge [ source 10 target 9 label \"-\" ]"
" edge [ source 11 target 9 label \"=\" ]"
" edge [ source 12 target 9 label \"-\" ]"
" edge [ source 13 target 12 label \"-\" ]"
" edge [ source 14 target 13 label \"-\" ]"
" edge [ source 15 target 14 label \"-\" ]"
" edge [ source 16 target 15 label \"-\" ]"
" edge [ source 17 target 16 label \"-\" ]"
" edge [ source 18 target 17 label \"-\" ]"
" edge [ source 19 target 17 label \"-\" ]"
" edge [ source 19 target 14 label \"-\" ]"
" edge [ source 20 target 19 label \"-\" ]"
" edge [ source 21 target 16 label \"-\" ]"
" edge [ source 22 target 21 label \":\" ]"
" edge [ source 23 target 22 label \":\" ]"
" edge [ source 24 target 23 label \":\" ]"
" edge [ source 25 target 24 label \":\" ]"
" edge [ source 26 target 25 label \"-\" ]"
" edge [ source 27 target 25 label \":\" ]"
" edge [ source 28 target 27 label \":\" ]"
" edge [ source 29 target 28 label \":\" ]"
" edge [ source 30 target 29 label \":\" ]"
" edge [ source 30 target 24 label \":\" ]"
" edge [ source 30 target 21 label \":\" ]"
" edge [ source 31 target 2 label \"-\" ]"
" edge [ source 32 target 31 label \"-\" ]"
" edge [ source 33 target 31 label \"-\" ]"
" edge [ source 33 target 0 label \"-\" ]"
" edge [ source 34 target 33 label \"-\" ]"
" edge [ source 2 target 35 label \"-\" ]"
" edge [ source 3 target 36 label \"-\" ]"
" edge [ source 3 target 37 label \"-\" ]"
" edge [ source 6 target 38 label \"-\" ]"
" edge [ source 10 target 39 label \"-\" ]"
" edge [ source 13 target 40 label \"-\" ]"
" edge [ source 13 target 41 label \"-\" ]"
" edge [ source 14 target 42 label \"-\" ]"
" edge [ source 16 target 43 label \"-\" ]"
" edge [ source 17 target 44 label \"-\" ]"
" edge [ source 18 target 45 label \"-\" ]"
" edge [ source 19 target 46 label \"-\" ]"
" edge [ source 20 target 47 label \"-\" ]"
" edge [ source 22 target 48 label \"-\" ]"
" edge [ source 26 target 49 label \"-\" ]"
" edge [ source 26 target 50 label \"-\" ]"
" edge [ source 28 target 51 label \"-\" ]"
" edge [ source 31 target 52 label \"-\" ]"
" edge [ source 32 target 53 label \"-\" ]"
" edge [ source 33 target 54 label \"-\" ]"
" edge [ source 34 target 55 label \"-\" ]"
"]"
	,
	"" // marks the end of the array
};



const std::string SMILES[] = {
	"[{CONH2}]c1ccc[n+](c1)[{Ribo-ADP}]",
	"[{CONH2}]C1[CH2]C=CN(C=1)[{Ribo-ADP}]",
	""
};


int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"    ggl::chem::MoleculeGroups  \n"
				<<"==============================================\n" 
				<<std::endl;

	  // container to fill from GML
	GroupMap groups;

	  // parse groups from GML
	for (size_t i = 0; !(gml[i].empty()); i++) {
		std::cout	<<"\n######### next GML #########\n"
					<<gml[i]
					<<"\n############################\n"
					<<std::endl;

		std::cout <<"\n-->MoleculeComponent_GML_grammar::parse( GML )" <<std::endl;

		  // do parsing
		std::pair<ggl::chem::MoleculeComponent, int>
				ret = MoleculeComponent_GMLparser::parseGML( gml[i] );

		if (ret.second < 0 )
		{
			  // print parsed group
			std::cout <<"\n parsed MoleculeComponent : '" <<ret.first.description <<"'" <<std::endl;
			std::cout <<Molecule_Graph(ret.first.pattern);
			std::cout <<"\n" <<std::endl;

			  // store group
			groups[ret.first.description] = ret.first;

		} else {
			std::cout <<"\n PARSING ERROR : at input position " <<ret.second
					<<" = '" <<gml[i].at(ret.second) <<"' "
					<<" within \n'"
					<<gml[i].substr(std::max(0,((int)ret.second)-10),20)
					<<"'\n"
					<<std::endl;
		}
	}

	  // parse SMILES and do tests
	for (size_t i = 0; !(SMILES[i].empty()); i++) {
		  // do parsing
		std::pair<ggl::chem::Molecule, int>
				mol = SMILESparser::parseSMILES( SMILES[i], groups );
		  // exception handling
		if (mol.second < 0 )
		{
			  // print parsed group
			std::cout <<"\n parsed Molecule SMILES : " <<SMILES[i] <<std::endl;

		} else {
			std::cout <<"\n PARSING ERROR : at input position " <<mol.second
					<<" = '" <<SMILES[i].at(mol.second) <<"' "
					<<" within \n'"
					<<SMILES[i].substr(std::max(0,((int)mol.second)-10),20)
					<<"'\n"
					<<std::endl;
			continue;
		}

		  // print without alteration
		std::cout <<"  new SMILES = " <<SMILESwriter::getSMILES(mol.first, true) <<std::endl;
		  // do group compression
		MoleculeUtil::compressGroups(mol.first, groups);
		std::cout <<" compression = " <<SMILESwriter::getSMILES(mol.first, groups, true) <<std::endl;

	}

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<std::endl;
	return 0;
}
