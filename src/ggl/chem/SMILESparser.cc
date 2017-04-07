

#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILES_grammar.hh"


#include <set>
#include <map>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "sgm/SubGraph.hh"
#include "sgm/RP_Hanser96.hh"

namespace ggl {
  namespace chem {
  
//##############################################################################
		


	std::pair< Molecule, int >
	SMILESparser
    ::parseSMILES( const std::string & SMILES_string )
	throw (std::invalid_argument)
    {
		  // forward call
		std::pair< Molecule, int > ret = SMILES_grammar::parseSMILES( SMILES_string );
		  // prune parsing artifacts if needed
		if ( ret.second == -1 ) {
			pruneAromaticNonRingBonds( ret.first );
		}
		return ret;
    }

//##############################################################################



	std::pair< Molecule, int >
	SMILESparser
    ::parseSMILES( const std::string & SMILES_string, const GroupMap& groups )
	throw (std::invalid_argument)
    {
		  // forward call
		std::pair< Molecule, int > ret = SMILES_grammar::parseSMILES( SMILES_string, groups );
		  // prune parsing artifacts if needed
		if ( ret.second == -1 ) {
			pruneAromaticNonRingBonds( ret.first );
		}
		return ret;
    }
	
//##############################################################################

    std::vector< Molecule * >
	SMILESparser::
    parseMultiSMILES( const std::string & SMILES_string )
    	throw (std::invalid_argument)
	{
    	GroupMap dummyMap;
    	return parseMultiSMILES( SMILES_string, dummyMap );
	}

//##############################################################################

    std::vector< Molecule * >
	SMILESparser::
    parseMultiSMILES( const std::string & SMILES_string, const GroupMap& groups )
    	throw (std::invalid_argument)
	{
    	if (SMILES_string.empty()) {
    		throw std::invalid_argument("empty multi SMILES string given for parsing");
    	}

		std::vector< Molecule* > mols;
		size_t firstPos = 0, endPos = 0;
		do {
			  // get next separator
			endPos = SMILES_string.find('.', firstPos);
			if (endPos == std::string::npos) {
				endPos = SMILES_string.size();
			}
			  // parse current SMILES
			std::pair< Molecule, int > mol = parseSMILES( SMILES_string.substr( firstPos, endPos - firstPos ), groups );

			if (mol.second == -1) {
				  // store molecule
				mols.push_back( new Molecule(mol.first) );
			} else {
				  // report error
				std::stringstream err;
				err <<"parsing error of SMILES '" <<SMILES_string <<"'"
						<<" at position " <<mol.second;
				throw std::invalid_argument(err.str());
			}

			firstPos = endPos + 1;
		} while ( endPos < SMILES_string.size() );

		  // return parsing result
		return mols;
	}

//##############################################################################

	ChemRule::CoreGraph
	SMILESparser::
	parseReactionSMILES( const std::string & SMILES_string, const bool pruneClassID )
	throw (std::invalid_argument)
	{
		 //////////////////  DECOMPOSE SMILES  ////////////////////

		  // get educts SMILES
		size_t firstPos = SMILES_string.find_first_not_of(" \t\n");
		if (firstPos == std::string::npos || SMILES_string.size() < 4) {
			throw std::invalid_argument("No reaction SMILES given");
		}
		size_t endPos = SMILES_string.find('>',firstPos);
		if (endPos == std::string::npos) {
			throw std::invalid_argument("Reaction SMILES misses first '>' separator");
		}
		const std::string eductSMILES = SMILES_string.substr(firstPos, endPos-firstPos);
		if (eductSMILES.empty()) {
			throw std::invalid_argument("Reaction SMILES has no educts");
		}

		firstPos = endPos + 1;
		endPos = SMILES_string.find('>', firstPos);
		if (endPos == std::string::npos) {
			throw std::invalid_argument("Reaction SMILES misses second '>' separator");
		}
		// ignore intermediate information

		  // get products
		firstPos = endPos + 1;
		if (firstPos >= SMILES_string.size()) {
			throw std::invalid_argument("Reaction SMILES has no products");
		}
		endPos = SMILES_string.find_first_of(" \t\n", firstPos);
		if (endPos == std::string::npos) {
			endPos = SMILES_string.size();
		}
		const std::string productSMILES = SMILES_string.substr(firstPos, endPos-firstPos);

		//////////////////  PARSE SMILES  ////////////////////

		Molecule educt, product;
		{
			std::vector< Molecule* > educts = parseMultiSMILES( eductSMILES );
			std::vector< Molecule* > products = parseMultiSMILES( productSMILES );

			// generate single educt / product graph covering all molecules

			MoleculeUtil::copy( Molecule_Graph_V( reinterpret_cast<std::vector< const Molecule* > &>(educts) ), educt );
			MoleculeUtil::copy( Molecule_Graph_V( reinterpret_cast<std::vector< const Molecule* > &>(products) ), product );

			// clear memory
			for (size_t i=0; i<educts.size(); ++i) {
				delete educts[i];
			} educts.clear();
			for (size_t i=0; i<products.size(); ++i) {
				delete products[i];
			} products.clear();
		}

		  // check for equal atom numbers
		if (boost::num_vertices(educt) != boost::num_vertices(product)) {
			throw std::invalid_argument("Reaction SMILES shows different numbers of atoms in educts and products");
		}


		//////////////////  get atom mapping  ////////////////////

		// vertex and edge iterator types
		boost::graph_traits<Molecule>::vertex_iterator     vi, vi_end;

		std::map<int,size_t> class2ruleID;
		size_t ruleID = 0;

		  // get atom mapping indices of educts
		std::vector<size_t> eductNewProtons(boost::num_vertices(educt), 0);
		std::vector<size_t> rule2educt(boost::num_vertices(educt), 0);
		std::vector<size_t> educt2rule(boost::num_vertices(educt), 0);
		boost::property_map<Molecule, PropNodeLabel>::type eNodeName = boost::get(PropNodeLabel(), educt);
		boost::property_map<Molecule, PropEdgeLabel>::type eEdgeName = boost::get(PropEdgeLabel(), educt);
		boost::property_map<Molecule, PropNodeIndex>::const_type eNodeID = boost::get(PropNodeIndex(), educt);
		for(boost::tie(vi, vi_end)=boost::vertices(educt); vi!=vi_end; ++vi) {
			const std::string label = eNodeName[*vi];
			int curID = MoleculeUtil::getClass( label );
			if (curID <= 0) {
				  // handle special case of adjacent protons
				if ( MoleculeUtil::getAtom(label) == "H" && boost::out_degree(*vi,educt) > 0) {
					  // get class ID from adjacent atom
					int neighID = MoleculeUtil::getClass( eNodeName[ boost::target(
																		*(boost::out_edges( *vi, educt ).first) /* first of one outgoing edge */
																		, educt) ] );
					  // check if ID
					if (neighID <= 0) {
						throw std::invalid_argument("Reaction SMILES of educts does not provide mapping ID in class for each non-proton atom");;
					}
					  // increase proton counter for this node
					++eductNewProtons[neighID];
					assert(eductNewProtons.at(neighID) < 10); // check needed for the following ID generation routine
					  // get the digit shift
					int digitShift = 10;
					while( neighID / digitShift > 0 ) {
						digitShift *= 10;
					}
					  // generate class ID for proton atom
					curID = ((int)rule2educt.size()) * digitShift * 10
							+ neighID * digitShift
							+ eductNewProtons.at(neighID);
				} else {
					throw std::invalid_argument("Reaction SMILES of educts does not provide mapping ID in class for each atom (except connected protons)");
				}
			} else if (class2ruleID.find(curID) != class2ruleID.end()) {
				std::stringstream err;
				err <<"Reaction SMILES : class ID '" <<curID <<"' is not unique in educts";
				throw std::invalid_argument(err.str());
			}
			  // store mapping
			class2ruleID[curID] = ruleID;
			++ruleID;
			rule2educt[class2ruleID[curID]] = eNodeID[*vi];
			educt2rule[eNodeID[*vi]] = class2ruleID[curID];
			  // insert newly generated IDs if not to be pruned
			if (!pruneClassID && curID > (int)rule2educt.size()) {
				const std::string atom = MoleculeUtil::getAtom(label);
				const int protons = MoleculeUtil::getProtons(label);
				const int charge = MoleculeUtil::getCharge(label);
				eNodeName[*vi] = MoleculeUtil::getComplexAtomLabel( atom, protons, charge, curID );

			}
		}
		if (pruneClassID) {
			for(boost::tie(vi, vi_end)=boost::vertices(educt); vi!=vi_end; ++vi) {
			  // prune class information from atom label
			  // TODO : maybe better to be replaced by MoleculeUtil functions
				eNodeName[*vi] = eNodeName[*vi].substr(0,eNodeName[*vi].find(':'));
//				const std::string atom = MoleculeUtil::getAtom(label);
//				const int protons = MoleculeUtil::getProtons(label);
//				const int charge = MoleculeUtil::getCharge(label);
//				eNodeName[*vi] = MoleculeUtil::getComplexAtomLabel( atom, protons, charge );

			}
		}

		  // get atom mapping indices of products
		std::set<int> classIDs;
		std::vector<size_t> productNewProtons(boost::num_vertices(educt), 0);
		std::vector<size_t> rule2product(boost::num_vertices(product), 0);
		std::vector<size_t> product2rule(boost::num_vertices(product), 0);
		boost::property_map<Molecule, PropNodeLabel>::type pNodeName = boost::get(PropNodeLabel(), product);
		boost::property_map<Molecule, PropEdgeLabel>::type pEdgeName = boost::get(PropEdgeLabel(), product);
		boost::property_map<Molecule, PropNodeIndex>::const_type pNodeID = boost::get(PropNodeIndex(), product);
		for(boost::tie(vi, vi_end)=boost::vertices(product); vi!=vi_end; ++vi) {
			const std::string label = pNodeName[*vi];
			int curID = MoleculeUtil::getClass( label );
			if (curID <= 0) {
				  // handle special case of adjacent protons
				if ( MoleculeUtil::getAtom(label) == "H" && boost::out_degree(*vi,product) > 0) {
					  // get class ID from adjacent atom
					int neighID = MoleculeUtil::getClass( pNodeName[ boost::target(
																		*(boost::out_edges( *vi, product ).first) /* first of one outgoing edge */
																		, product) ] );
					  // check if ID
					if (neighID <= 0) {
						throw std::invalid_argument("Reaction SMILES of products does not provide mapping ID in class for each non-proton atom");;
					}
					  // increase proton counter for this node
					++productNewProtons[neighID];
					assert(productNewProtons.at(neighID) < 10); // check needed for the following ID generation routine
					  // get the digit shift
					int digitShift = 10;
					while( neighID / digitShift > 0 ) {
						digitShift *= 10;
					}
					  // generate class ID for proton atom
					curID = ((int)rule2product.size()) * digitShift * 10
							+ neighID * digitShift
							+ productNewProtons.at(neighID);
					  // check if class ID is already known from educts
					if (class2ruleID.find(curID) != class2ruleID.end()) {
						throw std::invalid_argument("Reaction SMILES of products encodes an unbalanced proton without class ID");
					}
				} else {
					throw std::invalid_argument("Reaction SMILES of products does not provide mapping ID in class for each atom (except connected protons)");
				}
			} if (classIDs.find(curID) != classIDs.end()) {
				std::stringstream err;
				err <<"Reaction SMILES : class ID '" <<curID <<"' is not unique in products";
				throw std::invalid_argument(err.str());
			} if (class2ruleID.find(curID) == class2ruleID.end()) {
				std::stringstream err;
				err <<"Reaction SMILES : class ID '" <<curID <<"' of product atom is not present in educts";
				throw std::invalid_argument(err.str());
			}
			  // store mapping
			classIDs.insert(curID);
			rule2product[class2ruleID[curID]] = pNodeID[*vi];
			product2rule[pNodeID[*vi]] = class2ruleID[curID];
			  // insert newly generated IDs
			if (!pruneClassID && curID > (int)rule2educt.size()) {
				const std::string atom = MoleculeUtil::getAtom(label);
				const int protons = MoleculeUtil::getProtons(label);
				const int charge = MoleculeUtil::getCharge(label);
				pNodeName[*vi] = MoleculeUtil::getComplexAtomLabel( atom, protons, charge, curID );
			}
		}
		if (pruneClassID) {
			for(boost::tie(vi, vi_end)=boost::vertices(product); vi!=vi_end; ++vi) {
				  // prune class information from atom label
				  // TODO : maybe better to be replaced by MoleculeUtil functions
				pNodeName[*vi] = pNodeName[*vi].substr(0,pNodeName[*vi].find(':'));
			}
		}

		//////////////////  create rule core graph  ////////////////////

		ChemRule::CoreGraph rule;

		boost::property_map<ChemRule::CoreGraph, ChemRule::NodeContextProperty>::type
			nodeContext = boost::get( ChemRule::NodeContextProperty(), rule );
		boost::property_map<ChemRule::CoreGraph, ChemRule::NodeLabelProperty>::type
			nodeLabel = boost::get( ChemRule::NodeLabelProperty(), rule );
		boost::property_map<ChemRule::CoreGraph, ChemRule::NodeRightLabelProperty>::type
			nodeRightLabel = boost::get( ChemRule::NodeRightLabelProperty(), rule );
		boost::property_map<ChemRule::CoreGraph, ChemRule::EdgeContextProperty>::type
			edgeContext = boost::get( ChemRule::EdgeContextProperty(), rule );
		boost::property_map<ChemRule::CoreGraph, ChemRule::EdgeLabelProperty>::type
			edgeLabel = boost::get( ChemRule::EdgeLabelProperty(), rule );

		  // add nodes
		boost::graph_traits<ChemRule::CoreGraph>::vertex_descriptor rNode;
		for(size_t id = 0; id < rule2educt.size(); id++) {
			  // create new node
			rNode = boost::add_vertex( rule );
			  // set label
			nodeLabel[ rNode ] = eNodeName[ boost::vertex( rule2educt.at(id), educt ) ];
			  // check for label change
			if (nodeLabel[ rNode ] != pNodeName[ boost::vertex( rule2product.at(id), product ) ]) {
				  // set right label
				nodeRightLabel[ rNode ] = pNodeName[ boost::vertex( rule2product.at(id), product ) ];
				  // set label change context
				nodeContext[ rNode ] = ggl::Rule::RULE_LABEL_CHANGE;
			} else {
				  // no label change -> is a context node
				nodeContext[ rNode ] = ggl::Rule::RULE_CONTEXT;
			}
		}
		  // add left side and context edges
		boost::graph_traits<Molecule>::edge_iterator edge, edgeEnd;
		boost::graph_traits<ChemRule::CoreGraph>::edge_descriptor rEdge;
		for (boost::tie(edge,edgeEnd) = boost::edges( educt ); edge != edgeEnd; ++edge) {
			  // create new edge
			const size_t from = educt2rule.at( eNodeID[boost::source(*edge, educt)] );
			const size_t to   = educt2rule.at( eNodeID[boost::target(*edge, educt)] );
			rEdge = boost::add_edge( boost::vertex( from, rule)
									, boost::vertex( to, rule)
									, rule).first;
			  // set edge label
			edgeLabel[ rEdge ] = eEdgeName[ *edge ];
			  // check if left label differs from right one
			if (!boost::edge( boost::vertex(rule2product.at(from),product)
							, boost::vertex(rule2product.at(to),product)
							, product).second
				|| pEdgeName[boost::edge( boost::vertex(rule2product.at(from),product)
					, boost::vertex(rule2product.at(to),product)
					, product).first] != edgeLabel[ rEdge ] )
			{
				  // this is a true left side edge
				edgeContext[ rEdge ]  = ggl::Rule::RULE_LEFT_SIDE;
			} else {
				  // this is a context edge (same label educt/product)
				edgeContext[ rEdge ]  = ggl::Rule::RULE_CONTEXT;
			}
		}
		  // add missing right side edges
		for (boost::tie(edge,edgeEnd) = boost::edges( product ); edge != edgeEnd; ++edge) {
			  // create new edge
			const size_t from = product2rule.at( pNodeID[boost::source(*edge, product)] );
			const size_t to   = product2rule.at( pNodeID[boost::target(*edge, product)] );

			  // check if edge is not present or only for left side
			if ( !boost::edge(boost::vertex( from, rule)
							, boost::vertex( to, rule)
							, rule).second
				|| edgeContext[ boost::edge(boost::vertex( from, rule)
							, boost::vertex( to, rule)
							, rule).first ]  == ggl::Rule::RULE_LEFT_SIDE)
			{

				rEdge = boost::add_edge( boost::vertex( from, rule)
										, boost::vertex( to, rule)
										, rule).first;
				  // set edge label
				edgeLabel[ rEdge ] = pEdgeName[ *edge ];
				  // this is a true right side edge
				edgeContext[ rEdge ]  = ggl::Rule::RULE_RIGHT_SIDE;
			}
		}


		////////////////////  RETURN FINAL RULE GRAPH  ////////////////////

		return rule;
	}

//##############################################################################

	void
	SMILESparser::
	pruneAromaticNonRingBonds( Molecule& mol )
	{
		boost::property_map<	Molecule, PropNodeIndex >::const_type
			nodeIndex = boost::get( PropNodeIndex(), mol );
		boost::property_map<	Molecule, PropEdgeLabel >::type
			edgeLabel = boost::get( PropEdgeLabel(), mol );

		  // collect all aromatic bonds
		typedef sgm::RP_Hanser96::BondSet::value_type BOND;
		sgm::RP_Hanser96::BondSet aromBonds;
		Molecule::edge_iterator eIt, eItEnd;
		for (boost::tie(eIt,eItEnd) = boost::edges(mol); eIt != eItEnd; ++eIt) {
			  // check edge label if aromatic edge
			assert( MoleculeUtil::getBondData( edgeLabel[*eIt] ) != NULL);
			if ( MoleculeUtil::getBondData( edgeLabel[*eIt] )->isAromatic != 0) {
				  // store bond
				if( nodeIndex[boost::source(*eIt,mol)] < nodeIndex[boost::target(*eIt,mol)]) {
					aromBonds.insert( BOND( nodeIndex[boost::source(*eIt,mol)], nodeIndex[boost::target(*eIt,mol)]));
				} else {
					aromBonds.insert( BOND( nodeIndex[boost::target(*eIt,mol)], nodeIndex[boost::source(*eIt,mol)]));
				}
			}
		}

		  // check if any aromatic bonds
		if (aromBonds.empty()) {
			return;
		}

		  // ring bonds enumeration
		sgm::RP_Hanser96::BondSet ringBonds;
		sgm::RP_Hanser96 ringPerception;
		ringPerception.findRingBonds( Molecule_Graph(mol), ringBonds, 20 );

		  // remove all bonds part of ring
		std::vector<BOND> bondsToCorrect;
		std::set_difference( aromBonds.begin(), aromBonds.end()
							, ringBonds.begin(), ringBonds.end()
							, std::back_inserter(bondsToCorrect) );

		  // relabel remaining non-ring bonds to single bond
		Molecule::edge_descriptor curEdge;
		for (std::vector<BOND>::const_iterator b = bondsToCorrect.begin(); b != bondsToCorrect.end(); ++b) {
			  // access edge to rename
			curEdge = boost::edge( boost::vertex(b->first,mol)
								, boost::vertex(b->second,mol)
								, mol ).first;
			  // rename edge
			assert(edgeLabel[curEdge] == ":");
			edgeLabel[curEdge] = std::string("-");
		}

	}

//##############################################################################

  }
}
