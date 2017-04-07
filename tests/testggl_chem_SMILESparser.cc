
#include <iostream>
#include <iomanip>

#include "ggl/Graph.hh"
#include <sgm/Graph_boost.hh>

#include "ggl/Rule_GML_writer.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"

#include "dataTargetGraph_1.icc"

//#include "utilPrintGraph_Interface.icc"
#include "utilPrintRule.icc"


int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"         ggl::chem::SMILESparser  \n"
				<<"==============================================\n" 
				<<std::endl;
	
	std::vector<std::string> singleSMILES;
	singleSMILES.push_back( "COc1cccc(c1)C#N" );
	singleSMILES.push_back( "COc(cc1)ccc1C#N" );
	singleSMILES.push_back( "c1ccccc1c2ccccc2" );

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
				print(ggl::chem::Molecule_Graph( toFill.first ));
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
				<<"==============================================\n"
				<<std::endl;


	std::vector<std::string> multiSMILES;
	multiSMILES.push_back( "CC.O.[CH4]" );
	
	for (size_t i=0; i<multiSMILES.size(); i++) {
		
		std::cout <<"\n parse multi SMILES string " <<i <<" :\n" <<std::endl;
		
		std::cout <<" SMILES = '" <<multiSMILES[i] <<"'\n" <<std::endl;
	
		

		try {
			// parse smiles string into molecules
			std::vector< ggl::chem::Molecule*> toFill = ggl::chem::SMILESparser::parseMultiSMILES(multiSMILES.at(i));

			std::cout	<<" ==> multiSMILES Parser succeeded\n"
						<<" ==> resulting graph :\n"
						<<std::endl;
			  // print graph to stream
			print(ggl::chem::Molecule_Graph_V( reinterpret_cast<std::vector< const ggl::chem::Molecule* > &>(toFill) ));
			std::cout <<std::endl;

			// clear memory
			for (size_t m=0; m<toFill.size(); ++m) {
				delete toFill[m];
			}
			toFill.clear();

		} catch (std::exception &ex) {
			std::cout <<" EXCEPTION RAISED for parsing multiSMILES '"
						<< multiSMILES[i]
						<<"' :\n"
						<< ex.what();
		}
	
	}

	std::cout	<<"\n"
				<<"==============================================\n"
				<<std::endl;


	std::vector<std::string> reactionSMILES;
	reactionSMILES.push_back( "[C:13][C-2:2].[O:3]>>[O:3]=[C:2][C:13]" );

	for (size_t i=0; i<reactionSMILES.size(); i++) {

		std::cout <<"\n parse reaction SMILES string " <<i <<" :\n" <<std::endl;

		std::cout <<" SMILES = '" <<reactionSMILES[i] <<"'\n" <<std::endl;

		try {
			// parse smiles string into molecules
			ggl::chem::ChemRule::CoreGraph toFill = ggl::chem::SMILESparser::parseReactionSMILES(reactionSMILES.at(i));

			std::cout	<<" ==> reactionSMILES Parser succeeded\n"
						<<" ==> create rule :\n"
						<<std::endl;

			ggl::chem::ChemRule rule(toFill, "parsed rule");

			ggl::Rule_GML_writer::write(std::cout, rule, true );

			  // print rule to stream
			printRule(rule);
			std::cout <<std::endl;

		} catch (std::exception &ex) {
			std::cout <<" EXCEPTION RAISED for parsing reactionSMILES '"
						<< reactionSMILES[i]
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
