
#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILESwriter.hh"
#include "ggl/chem/MoleculeUtil.hh"

#include <sgm/Graph_boost.hh>


#include "utilPrintGraph_Interface.icc"

using namespace ggl::chem;

#include <cassert>
#include <algorithm>

/**
 * creates a node shuffled copy of the given molecule according to the index map.
 */
Molecule getShuffle( const Molecule & mol, const std::vector<size_t> indexMap ) {
	assert(boost::num_vertices(mol) == indexMap.size());

	boost::property_map< Molecule, PropNodeLabel>::const_type
		mNodeLabel = (boost::get( PropNodeLabel(), mol ));
	boost::property_map< Molecule, PropEdgeLabel>::const_type
		mEdgeLabel = (boost::get( PropEdgeLabel(), mol ));
	boost::property_map< Molecule, PropNodeIndex>::const_type
		mNodeIndex = (boost::get( PropNodeIndex(), mol ));

	Molecule ret;
	boost::property_map< Molecule, PropNodeLabel>::type
		rNodeLabel = (boost::get( PropNodeLabel(), ret ));
	boost::property_map< Molecule, PropEdgeLabel>::type
		rEdgeLabel = (boost::get( PropEdgeLabel(), ret ));
	boost::property_map< Molecule, PropNodeIndex>::const_type
		rNodeIndex = (boost::get( PropNodeIndex(), ret ));

	  // get reverse mapping
	const std::vector<size_t>& mol2ret = indexMap;
	std::vector<size_t> ret2mol(mol2ret.size(),0);
	for (size_t i=0; i<mol2ret.size(); ++i) {
		ret2mol[mol2ret.at(i)] = i;
	}
	  // add nodes and set mapped label
	for (size_t i=0; i<ret2mol.size(); ++i) {
		rNodeLabel[boost::add_vertex(ret)] = mNodeLabel[boost::vertex(ret2mol.at(i),mol)];
	}
	  // add edges and set labels
	Molecule::edge_iterator ei, e_end;
	for (boost::tie(ei,e_end) = boost::edges(mol); ei!=e_end; ++ei) {
		 // add edge
		Molecule::edge_descriptor edge =
			boost::add_edge(
			mol2ret.at(mNodeIndex[boost::source(*ei,mol)])
			, mol2ret.at(mNodeIndex[boost::target(*ei,mol)])
			, ret ).first;
		 // set label
		rEdgeLabel[edge] = mEdgeLabel[*ei];
	}
	  // final molecule returned
	return ret;
}


int main() {

	std::vector<std::string> SMILES;
	
	  // basics
	SMILES.push_back("[H+]");
	SMILES.push_back("HH");
	SMILES.push_back("[HH]");
	SMILES.push_back("[H][H]");
	SMILES.push_back("O");
	SMILES.push_back("[OH-]");
	SMILES.push_back("[O-]");
	SMILES.push_back("[O+]");
	SMILES.push_back("[OH3+]");


	SMILES.push_back( "CC(C)C(=O)O" );
	SMILES.push_back( "n1ccccc1" );
	
	SMILES.push_back( "[nH]1cccc1" );
	SMILES.push_back( "Hn1cccc1" );
	
	SMILES.push_back( "CC(C)(N)Cc1ccccc1" );

	SMILES.push_back( "C(OP([O-])(=O)[O-])C1(OC(O)CC1(O))" );

	// equal
	SMILES.push_back( "Cn1cnc2n(C)c(=O)n(C)c(=O)c12" );
	SMILES.push_back( "Cn1cnc2c1c(=O)n(C)c(=O)n2C" );
	
	// equal
	SMILES.push_back( "CCN(CC)C(=O)" );
	SMILES.push_back( "CCN(C=O)CC" );
	
	// equal 
	SMILES.push_back( "COc1cccc(C#N)c1" );
	SMILES.push_back( "COc(c1)cccc1C#N" );
	SMILES.push_back( "COc1cccc(CN)c1" );
	
	// equal
	SMILES.push_back( "COc(cc1)ccc1C#N" );
	SMILES.push_back( "COc1ccc(CN)cc1" );
	
	// equal
	SMILES.push_back( "c1ccccc1-c2ccccc2" );
	SMILES.push_back( "c1ccc(cc1)-c2ccccc2" );

	// with coordination bonds
	SMILES.push_back( "C^C" );
	SMILES.push_back( "[Fe+2](^O)^O" );
	SMILES.push_back( "c1cc(O2)c([O-]^[Fe+2]^2)cc1" );
	
	// equal
	SMILES.push_back( "[C][C:1]" );
	SMILES.push_back( "[C:1][C]" );

	// equal
	SMILES.push_back( "C1C[C:1][C:2]=[C:3][C:4]1" );
	SMILES.push_back( "C1C[C:4][C:3]=[C:2][C:1]1" );

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
						<< "stopped at: \"" << SMILES[i][result.second] << "\"" <<"\n"
						<< "\"" << SMILES[i] << "\"" <<"\n"
						<< std::setw(1+result.second) << " " << "^ error"
						<<"\n"
						<< std::endl;
			  // go to next
			continue;
		}
		

		std::cout <<" --> call newSMILES = SMILESwriter::getSMILES(..)  = "; std::cout.flush();
		std::string newSMILES = ggl::chem::SMILESwriter::getSMILES(result.first);
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;
		
		std::cout <<" --> call MoleculeUtil::fillProtons(..) "<<std::endl;
		ggl::chem::MoleculeUtil::fillProtons( result.first );
		std::cout <<"  > call newSMILES = SMILESwriter::getSMILES( filled, ignoreProtons=false ) = "; std::cout.flush();
		newSMILES = ggl::chem::SMILESwriter::getSMILES( result.first, false );
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;
		std::cout <<"  > call newSMILES = SMILESwriter::getSMILES( filled, ignoreProtons=true )  = "; std::cout.flush();
		newSMILES = ggl::chem::SMILESwriter::getSMILES( result.first, true );
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;

		std::cout <<" --> call MoleculeUtil::compressHnodes( compMol ) "<<std::endl;
		ggl::chem::MoleculeUtil::compressHnodes( result.first );
		std::cout <<"  > call newSMILES = SMILESwriter::getSMILES( compMol, ignoreProtons=false) = "; std::cout.flush();
		newSMILES = ggl::chem::SMILESwriter::getSMILES( result.first, false );
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;
		std::cout <<"  > call newSMILES = SMILESwriter::getSMILES( compMol, ignoreProtons=true ) = "; std::cout.flush();
		newSMILES = ggl::chem::SMILESwriter::getSMILES( result.first, true );
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;
		
		std::cout <<" --> call MoleculeUtil::removeProtons(..) "<<std::endl;
		ggl::chem::MoleculeUtil::removeProtons( result.first );
		std::cout <<"  > call newSMILES = SMILESwriter::getSMILES( remMol, ignoreProtons=false ) = "; std::cout.flush();
		newSMILES = ggl::chem::SMILESwriter::getSMILES( result.first, false );
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;
		std::cout <<"  > call newSMILES = SMILESwriter::getSMILES( remMol, ignoreProtons=true )  = "; std::cout.flush();
		newSMILES = ggl::chem::SMILESwriter::getSMILES( result.first, true );
		std::cout <<"'" <<newSMILES <<"'"<<std::endl;

	}

	SMILES.clear();
	SMILES.push_back( "[H][C]12CCC3=Cc4oncc4C[C]3(C)[C]1([H])CC[C]1(C)[C]2([H])CC[C]1(O)C#C" );



	for (size_t i=0; i<SMILES.size(); i++) {

		std::cout <<"\n\n PERMUTATING SMILES = '" <<SMILES[i] <<"'" <<std::endl;

//		std::cout <<" --> call SMILESparser::parse(..) "<<std::endl;
		std::pair<ggl::chem::Molecule,int> result
			= ggl::chem::SMILESparser::parseSMILES(SMILES[i]);
		if (result.second != -1) {
			std::cout	<< " ==> SMILES Parser failed !!!\n"
						<< "stopped at: \"" << SMILES[i][result.second] << "\"" <<"\n"
						<< "\"" << SMILES[i] << "\"" <<"\n"
						<< std::setw(1+result.second) << " " << "^ error"
						<<"\n"
						<< std::endl;
			  // go to next
			continue;
		}

		size_t numNonProtons = boost::num_vertices(result.first);
		ggl::chem::MoleculeUtil::fillProtons( result.first );

		std::set<std::string> permutations;

		std::vector<size_t> indexMap(boost::num_vertices(result.first),0);
		for(size_t i=0; i<indexMap.size(); ++i) { indexMap[i] = i; }
		size_t permNumber = 0;

		for (size_t first=0; first<numNonProtons; ++first) {
			std::vector<size_t> rotMap = indexMap;
			std::rotate_copy( indexMap.begin(), indexMap.begin()+first, indexMap.begin()+numNonProtons, rotMap.begin());
			Molecule permutation = getShuffle( result.first, rotMap );
//
//		do {
//			Molecule permutation = getShuffle( result.first, indexMap );
			std::string permSMILES = ggl::chem::SMILESwriter::getSMILES( permutation );
			if ( permutations.find(permSMILES) == permutations.end()) {
				std::cout <<"   permutation number "<<permNumber <<" : " << permSMILES <<std::endl;
//				std::cout <<Molecule_Graph(permutation)<<std::endl;
				permutations.insert(permSMILES);
			}
			++permNumber;
		}
//		while( std::next_permutation(indexMap.begin(),indexMap.begin()+numNonProtons) );
	}


	return 0;
}


