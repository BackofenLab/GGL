
#include <iostream>
#include <iomanip>

#include "ggl/Graph.hh"
#include <sgm/Graph_boost.hh>

#include "ggl/Rule_GML_writer.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/Molecule_Graph_noClass.hh"

#include "dataTargetGraph_1.icc"

//#include "utilPrintGraph_Interface.icc"
#include "utilPrintRule.icc"


int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"         ggl::chem::Graph_noClass\n"
				<<"==============================================\n" 
				<<std::endl;
	
	std::vector<std::string> singleSMILES;
	singleSMILES.push_back( "[CH3+:123]" );
	singleSMILES.push_back( "[H:1][C:2]#[N:3]" );

	for (size_t i=0; i<singleSMILES.size(); i++) {

		std::cout <<"\n parse single SMILES string " <<i <<" :\n" <<std::endl;

		std::cout <<" SMILES = '" <<singleSMILES[i] <<"'\n" <<std::endl;



		try {
			// parse smiles string into molecules
			std::pair< ggl::chem::Molecule, int> toFill = ggl::chem::SMILESparser::parseSMILES(singleSMILES.at(i));
			if (toFill.second == -1) {

				std::cout	<<" ==> singleSMILES Parser succeeded\n"
							<<" ==> resulting graph :\n"
							<<std::endl;
				  // print graph to stream
				ggl::chem::Molecule_Graph mol( toFill.first );
				print(mol);
				std::cout	<<" ==> without class :\n"
							<<std::endl;
				ggl::chem::Molecule_Graph_noClass molNoClass( mol );
				print(molNoClass);
				std::cout	<<" ==> original molecule :\n"
							<<std::endl;
				print(molNoClass.getWithFullAtomLabels());
				std::cout <<std::endl;
			} else {
				throw std::exception();
			}


		} catch (std::exception &ex) {
			std::cout <<" EXCEPTION RAISED for parsing singleSMILES '"
						<< singleSMILES[i]
						<<"' :\n"
						<< ex.what();
		}

	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
