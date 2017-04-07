
#include <iostream>
#include <iomanip>

#include "ggl/Graph.hh"
#include <sgm/Graph_boost.hh>

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILES_grammar.hh"

#include "dataTargetGraph_1.icc"

#include "utilPrintGraph_Interface.icc"

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"         ggl::chem::SMILES_grammar  \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	std::vector<std::string> SMILES;
	
	SMILES.push_back( "CC(C)C(=O)O" );
	SMILES.push_back( "n1ccccc1" );
	SMILES.push_back( "[nH]1cccc1" );
	SMILES.push_back( "[CH3][CH]=[CH][CH3]" );

	 // all the following 18 SMILES encode the same molecule
	SMILES.push_back( "C=12CC1C2" );
	SMILES.push_back( "C=12CC1C-2" );
	SMILES.push_back( "C=1-2CC1C2" );
	SMILES.push_back( "C=1-2CC1C-2" );
	SMILES.push_back( "C12CC=1C2" );
	SMILES.push_back( "C12CC=1C-2" );
	SMILES.push_back( "C1-2CC=1C2" );
	SMILES.push_back( "C1-2CC=1C-2" );
	SMILES.push_back( "C=1-2CC=1C-2" );
	SMILES.push_back( "C1=2CC2C1" );
	SMILES.push_back( "C-1=2CC2C1" );
	SMILES.push_back( "C1=2CC2C-1" );
	SMILES.push_back( "C-1=2CC2C-1" );
	SMILES.push_back( "C12CC=2C1" );
	SMILES.push_back( "C-12CC=2C1" );
	SMILES.push_back( "C12CC=2C-1" );
	SMILES.push_back( "C-12CC=2C-1" );
																																			SMILES.push_back( "C-1=2CC=2C-1" );
	SMILES.push_back( "O[As]" );
	SMILES.push_back( "[H][C]12CCC3=Cc4oncc4C[C]3(C)[C]1([H])CC[C]1(C)[C]2([H])CC[C]1(O)C#C" );
	SMILES.push_back( "CC(O)C(O)C1CNc2nc(N)[nH]c(=O)c2N1" );
	SMILES.push_back( "N#CSCc1ccccc1" );
	SMILES.push_back( "CN(C)c1ccc2nc3ccc(cc3[s+]c2c1)N(C)C" );
	SMILES.push_back( "C1=CC=[N](C=C1)[Ag++]([N]1=CC=CC=C1)([N]1=CC=CC=C1)[N]1=CC=CC=C1" );

	SMILES.push_back( "c1ccccc1c2ccccc2" );
	SMILES.push_back( "c1ccccc1-c2ccccc2" );
	SMILES.push_back( "c1ccc(cc1)-c2ccccc2" );

	SMILES.push_back( "c1cc[se]c1" );

	SMILES.push_back( "[c:13]1cc[se]c1C[Br:5]" );

	 // with coordination bonds
	SMILES.push_back( "C^C" );
	SMILES.push_back( "[Fe+2](^O)^O" );

	 // special ring closure case (with leading bond)
	SMILES.push_back( "C(C%10CC=%10)=C" );



	for (size_t i=0; i<SMILES.size(); i++) {
		
		std::cout <<"\n convert SMILES string " <<i <<" :\n" <<std::endl;
		
		std::cout <<" SMILES = '" <<SMILES[i] <<"'\n" <<std::endl;
	
		 // the graph to fill
		ggl::chem::Molecule toFill;
	
		 // create SMILES grammar parser
		NS_BOOSTSPIRIT::parse_info<> info;
		ggl::chem::SMILES_grammar smiles_grammar(toFill);
		

		try {
			// parse smiles string into graph
			info = NS_BOOSTSPIRIT::parse(SMILES[i].c_str(), smiles_grammar);

			// check if parsing was successful
			if (info.full) {
				std::cout	<<" ==> SMILES Parser succeeded\n"
							<<" ==> resulting graph :\n"
							<<std::endl;
				  // print graph to stream
				sgm::Graph_boost<ggl::Graph> toPrint(toFill);
				print(toPrint);
				std::cout <<std::endl;
			}
			else {
				  // give parsing error output
				std::cout	<< " ==> SMILES Parser failed !!!\n"
							<< "stoped at: \"" << info.stop << "\"" <<"\n"
							<< "\"" << SMILES[i] << "\"" <<"\n"
							<< std::setw(1+info.length) << " " << "^ error"
							<<"\n"
							<< std::endl;
			}

			std::cout <<" --> call SMILES_grammar::parse(..) "<<std::endl;
			std::pair<ggl::Graph,int> result
				= ggl::chem::SMILES_grammar::parseSMILES(SMILES[i]);
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
			}
		} catch (std::exception &ex) {
			std::cout <<" EXCEPTION RAISED for parsing SMILES '"
						<< SMILES[i]
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
