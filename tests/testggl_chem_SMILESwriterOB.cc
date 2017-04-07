
#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILESwriterOB.hh"

#include <sgm/Graph_boost.hh>


#include "utilPrintGraph_Interface.icc"

#include <cstring>

int main() {

	std::vector<std::string> SMILES;
	
	SMILES.push_back( "CC(C)C(=O)O" );
	SMILES.push_back( "N1CCCCC1" );
	SMILES.push_back( "C1NCCCC1" );
	SMILES.push_back( "C1CNCCC1" );
	SMILES.push_back( "C1CCNCC1" );
	SMILES.push_back( "C1CCCNC1" );
	SMILES.push_back( "C1CCCCN1" );
	SMILES.push_back( "[NH]1CCCC1" );
	
	
	for (size_t i=0; i<SMILES.size(); i++) {

		std::cout <<"\n\n PROCESSING SMILES = '" <<SMILES[i] <<"'\n" <<std::endl;

		std::cout <<" --> call SMILESparser::parse(..) "<<std::endl;
		std::pair<ggl::chem::Molecule,int> result 
			= ggl::chem::SMILESparser::parseSMILES(SMILES[i]);
		if (result.second == -1) {
			std::cout	<<" ==> SMILES Parser succeeded\n"
						<<" ==> resulting graph :\n"
						<<std::endl;
			  // print graph to stream
			sgm::Graph_boost<ggl::Graph> toPrint(result.first);
			print(toPrint);
			std::cout <<std::endl;
		} else {
			std::cout	<< " ==> SMILES Parser failed !!!\n"
						<< "stoped at: \"" << SMILES[i][result.second] << "\"" <<"\n"
						<< "\"" << SMILES[i] << "\"" <<"\n"
						<< std::setw(1+result.second) << " " << "^ error"
						<<"\n"
						<< std::endl;
			  // go to next
			continue;
		}
		
		std::cout <<" --> call newSMILES = SMILESwriterOB::getSMILES(..) "<<std::endl;
		std::string newSMILES = ggl::chem::SMILESwriterOB::getSMILES(result.first);
		std::cout <<" newSMILES = '" <<newSMILES <<"'"<<std::endl;

	}
	return 0;
}
