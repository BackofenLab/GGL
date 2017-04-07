

#include <iostream>
#include <iomanip>

#include "ggl/Graph.hh"
#include <sgm/Graph_boost.hh>

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/GS_SMILES.hh"

#include "dataTargetGraph_1.icc"

#include "utilPrintGraph_Interface.icc"

void runSuccessiveInsert( const std::vector<std::string>& SMILES) {
	
	typedef std::set< std::string > SMILES_container;
	typedef std::insert_iterator< SMILES_container > SMILES_inserter;
	
//	typedef std::vector< std::string > SMILES_container;
//	typedef std::back_insert_iterator< SMILES_container > SMILES_inserter;

//	typedef std::list< std::string > SMILES_container;
//	typedef std::front_insert_iterator< SMILES_container > SMILES_inserter;
	
	SMILES_container smilesStorage;
	
	SMILES_inserter insert(smilesStorage,smilesStorage.begin()); // for insert_iterator
//	SMILES_inserter insert(smilesStorage);	// for front_* or back_insert_iterator
	
	std::cout	<<"\n create GS_SMILES(INSERTER) gs"<<std::endl;
	ggl::chem::GS_SMILES< SMILES_inserter > gs(insert);
	
	for (size_t i=0; i<SMILES.size(); i++) {
		
		std::cout <<"\n convert SMILES string " <<i <<" :\n" <<std::endl;
		
		std::cout <<" SMILES = '" <<SMILES.at(i) <<"'\n" <<std::endl;
	
		std::cout <<" --> call SMILESparser::parse(..) "<<std::endl;
		std::pair<ggl::Graph,int> result 
			= ggl::chem::SMILESparser::parseSMILES(SMILES.at(i));
		if (result.second == -1) {
			std::cout	<<" ==> SMILES Parser succeeded\n"
						<<" ==> resulting graph :\n"
						<<std::endl;
			  // print graph to stream
			sgm::Graph_boost<ggl::Graph> toPrint(result.first);
			print(toPrint);
			std::cout <<std::endl;
			std::cout <<"\n ==> adding graph to gs"<<std::endl;
			gs.add(result.first);
			std::cout <<"\n ==> adding graph to gs"<<std::endl;
			gs.add(result.first);

		} else {
			std::cout	<< " ==> SMILES Parser failed !!!\n"
						<< "stoped at: \"" << SMILES.at(i)[result.second] << "\"" <<"\n"
						<< "\"" << SMILES.at(i) << "\"" <<"\n"
						<< std::setw(1+result.second) << " " << "^ error"
						<<"\n"
						<< std::endl;
		}
	}
	
	std::cout <<" ==> content of SMILES storage :" <<std::endl;
	for ( SMILES_container::const_iterator it = smilesStorage.begin();
			it != smilesStorage.end(); it++ )
	{
		std::cout	<<"    " <<*it <<std::endl;
	}

}

void runFullInsert( const std::vector<std::string>& SMILES) {
	
	typedef std::set< std::string > SMILES_container;
	typedef std::insert_iterator< SMILES_container > SMILES_inserter;
	
//	typedef std::vector< std::string > SMILES_container;
//	typedef std::back_insert_iterator< SMILES_container > SMILES_inserter;

//	typedef std::list< std::string > SMILES_container;
//	typedef std::front_insert_iterator< SMILES_container > SMILES_inserter;
	
	SMILES_container smilesStorage;
	
	SMILES_inserter insert(smilesStorage,smilesStorage.begin()); // for insert_iterator
//	SMILES_inserter insert(smilesStorage);	// for front_* or back_insert_iterator
	
	std::cout	<<"\n create GS_SMILES(INSERTER) gs"<<std::endl;
	ggl::chem::GS_SMILES< SMILES_inserter > gs(insert);
	
	 // the graph to fill
	ggl::chem::Molecule toFill;

	for (size_t i=0; i<1 && i<SMILES.size(); i++) {
		
		std::cout <<"\n convert SMILES string " <<i <<" :\n" <<std::endl;
		
		std::cout <<" SMILES = '" <<SMILES.at(i) <<"'\n" <<std::endl;
	
		 // create SMILES grammar parser
		std::pair< ggl::chem::Molecule, int > ret =
				ggl::chem::SMILESparser::parseSMILES( SMILES.at(i) );

		// check if parsing was successful
		if (ret.second < 0) {

			  // make copy for further processing
			ggl::chem::MoleculeUtil::copy( ret.first, toFill );

			std::cout	<<" ==> SMILES Parser succeeded\n"
						<<" ==> resulting graph :\n"
						<<std::endl;
			  // print graph to stream
			sgm::Graph_boost<ggl::Graph> toPrint(toFill);
			print(toPrint);
		}
		else {
			  // give parsing error output
			std::cout	<< " ==> SMILES Parser failed !!!\n"
						<< "stopped at: \"" << ret.second << "\"" <<"\n"
						<< "\"" << SMILES.at(i) << "\"" <<"\n"
						<< std::endl;
		}
	}
	
	boost::property_map<ggl::chem::Molecule, ggl::PropNodeLabel>::type MV_Property_Map_t = boost::get(ggl::PropNodeLabel(), toFill);
	boost::property_map<ggl::chem::Molecule, ggl::PropEdgeLabel>::type ME_Property_Map_t = boost::get(ggl::PropEdgeLabel(), toFill);
	
	ggl::chem::Molecule::vertex_descriptor c1 = boost::add_vertex(toFill);
	MV_Property_Map_t[c1] = "C";
	ggl::chem::Molecule::vertex_descriptor c2 = boost::add_vertex(toFill);
	MV_Property_Map_t[c2] = "C";
	ggl::chem::Molecule::vertex_descriptor o1 = boost::add_vertex(toFill);
	MV_Property_Map_t[o1] = "O";
	ggl::chem::Molecule::edge_descriptor e1 = boost::add_edge(c1, c2, toFill).first;
	ME_Property_Map_t[e1] = "=";
	ggl::chem::Molecule::edge_descriptor e2 = boost::add_edge(c2, o1, toFill).first;
	ME_Property_Map_t[e2] = "-";
	
	std::cout	<<"\n ==> final graph :\n"
				<<std::endl;
	  // print graph to stream
	sgm::Graph_boost<ggl::Graph> toPrint(toFill);
	print(toPrint);

	std::cout <<"\n ==> adding graph to gs"<<std::endl;
	gs.add(toFill);

	
	std::cout <<" ==> content of SMILES storage :" <<std::endl;
	for ( SMILES_container::const_iterator it = smilesStorage.begin();
			it != smilesStorage.end(); it++ )
	{
		std::cout	<<"    " <<*it <<std::endl;
	}

}

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"         ggl::chem::GS_SMILES  \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	std::vector<std::string> SMILES;
	
	SMILES.push_back( "CC(C)C(=O)O" );
	SMILES.push_back( "n1ccccc1" );
	SMILES.push_back( "[nH]1cccc1" );
	
	runSuccessiveInsert(SMILES);
	runFullInsert(SMILES);

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
