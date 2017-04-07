
#include <cassert>
#include <cmath>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <ostream>

#include <boost/graph/copy.hpp>

#include "ggl/chem/MoleculeUtil.hh"
#include "ggl/chem/AromaticityPerception.hh"


namespace ggl {
namespace chem {


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	bool
	AromaticityPerception::AdjacencyComp::
	operator() ( AromaticityPerception::AdjacencyData* e1, AromaticityPerception::AdjacencyData* e2)  {

		assert (e1 != NULL);
		assert (e2 != NULL);

		// check if (e1 smaller e2)

		 // check for equality
		bool bothSameOpenEdges = e1->openEdges == e2->openEdges;
		bool bothSameRemValence = e1->remValence == e2->remValence;

		 // less conditions
		bool e1NoOpenEdge = e1->openEdges == 0;
		bool e2NoOpenEdge = e2->openEdges == 0;
		bool e1LessOpenEdges = e1->openEdges < e2->openEdges;
		bool e1LessRemValence = e1->remValence < e2->remValence;
		bool e1LessAromEdges = e1->aromaticEdges < e2->aromaticEdges;
		bool e1OnlySingleEdges = e1->openEdges == e1->remValence;
		bool e2OnlySingleEdges = e2->openEdges == e2->remValence;

		////////////////// e1/e2 has no edges left

		  // no open edge check
		if (e1NoOpenEdge || e2NoOpenEdge) {
			  // e1 has to have no open edge to be smaller
			return	(!e1NoOpenEdge && e2NoOpenEdge);
		}

		////////////  both have at least 1 open edge

		  // single edges are smallest
		if (e1OnlySingleEdges && ! e2OnlySingleEdges) {
			return true;
		}

		  // single edges tie break via remaining valence
		if (e1OnlySingleEdges && e2OnlySingleEdges) {
			return e1LessRemValence || (bothSameRemValence && e1LessAromEdges);
		}

		  // final sorting via number of open edges
		return	e1LessOpenEdges
				||
				(bothSameOpenEdges && e1LessRemValence)
				||
				(bothSameOpenEdges && bothSameRemValence && e1LessAromEdges)
				;

	}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	void
	AromaticityPerception::
	correctAromaticity( Molecule & mol, const bool checkValence )
	 throw (std::runtime_error)
	{
		  // clear temporary data
		clearData();
		  // identify all rings
		findAllRings( mol );

		  // prune rings to remove fused representation
		pruneFusedRings( allRings );

		  // prune rings that cannot be aromatic
		pruneNonSingleDoubleBondRings( allRings, mol );

		  // identify all aromatic rings
		identifyAromaticEdges(mol);

		  // relabel according to identified aromatic rings
		relabelMolecule(mol, checkValence);
	}


////////////////////////////////////////////////////////////////////////////////



	void
	AromaticityPerception::
	pruneNonSingleDoubleBondRings( std::vector< RingDescriptor* > & rings
					, const Molecule & mol )
	{

		  // access to edge label
		boost::property_map<	Molecule , PropEdgeLabel >
			::const_type edgeLabel = boost::get( PropEdgeLabel(), mol );

		const MoleculeUtil::BondLabelData * bond = NULL;

		  // check each ring
		for (size_t i=0; i<rings.size(); ) {

			bool remove = false;
			  // check each bond if not single, double, or aromatic bond
			for (EdgeSet::const_iterator e = rings.at(i)->edges.begin(); !remove && e != rings.at(i)->edges.end(); ++e) {
				  // get bond data
				bond = MoleculeUtil::getBondData( edgeLabel[
										  boost::edge( boost::vertex(e->first,mol), boost::vertex(e->second,mol), mol).first
				                                                                          ] );
				assert(bond != NULL);
				remove = bond->valence > 2 || bond->valence < 1;
			}
			  // check if this ring has to be removed
			if (remove) {
				  // delete element
				delete rings[i];
				  // delete entry
				rings.erase(rings.begin()+i);
			} else {
				  // go to next element
				++i;
			}
		}

	}


////////////////////////////////////////////////////////////////////////////////



	void
	AromaticityPerception::
	pruneFusedRings( std::vector< RingDescriptor* > & rings )
	{

		  // sort rings by size
		RingSizeLess ringSizeLess;
		std::sort( rings.begin(), rings.end(), ringSizeLess );

//std::cerr <<"\nNEXT\n";
		  // prune rings
		size_t numOfprunedRings = 0;
		for (int i=rings.size()-1; i>=0; i--) {
//std::cerr <<"next ring = ";
//for (RingList::const_iterator x=rings.at(i)->ring.begin(); x!=rings.at(i)->ring.end(); ++x) {
//	std::cerr <<*x <<" ";
//}
			bool canBePruned = false;
			for( size_t c=0; !canBePruned && c < (size_t)i; ++c ) {
				  // skip already deleted rings
				if (rings.at(c) == NULL) {
					continue;
				}
				  // stop if ring is of equal size or larger
				if (rings.at(c)->edges.size() >= rings.at(i)->edges.size()) {
					break;
				}
				  // check for containment of this ring c within ring i

				  // get iterators on edge sets to compare
				EdgeSet::const_iterator largeEdge = rings.at(i)->edges.begin(), largeEnd = rings.at(i)->edges.end();
				EdgeSet::const_iterator curEdge = rings.at(c)->edges.begin(), curEnd = rings.at(c)->edges.end();
				  // sequential check
				size_t curOverlap = 0;
				  // check only until up to two differences have been found
				while(largeEdge != largeEnd && curEdge != curEnd) {
					  // check if overlap
					if (*largeEdge == *curEdge) {
						++curOverlap;
						++largeEdge;
						++curEdge;
					  // check which counter to increase
					} else if (*largeEdge < *curEdge) {
						++largeEdge;
					} else {
						++curEdge;
					}
				}
				  // check if smaller current ring is contained excluding one edge
				  // if only one edge is left : ring is contained
				canBePruned = (curOverlap+1 == rings.at(c)->edges.size());

			}
			  // prune ring if obsolete
			if (canBePruned) {
				  // remove ring descriptor at position i
				delete rings[i];
				  // set to NULL to remember deletion
				rings[i] = NULL;
				  // count pruning
				numOfprunedRings++;
//std::cerr <<" pruned";
			}
//std::cerr <<"\n";
		}

		  // shift remaining rings to the front
		if (numOfprunedRings>0) {
			size_t fillPos = 0;
			for( size_t i=0; i<rings.size(); ++i ) {
				if (rings.at(i) != NULL) {
					  // check if overwrite is needed
					if (fillPos < i) {
						rings[fillPos] = rings[i];
					}
					  // increase overwrite counter
					fillPos++;
				}
			}
			  // drop deleted elements
			rings.resize(rings.size()-numOfprunedRings);
		}

	}

////////////////////////////////////////////////////////////////////////////////



	void
	AromaticityPerception::
	relabelMolecule( Molecule & result, const bool checkValence )
	 throw (std::runtime_error)
	{
//		  // new molecule to be filled
//		Molecule result(mol);

		  // property maps for access
		boost::property_map<	Molecule , PropNodeLabel >
			::type resultNodeLabel = boost::get( PropNodeLabel(), result );

		boost::property_map<	Molecule , PropNodeIndex >
			::type resultNodeIndex = boost::get( PropNodeIndex(), result );

		boost::property_map<	Molecule , PropEdgeLabel >
			::type resultEdgeLabel = boost::get( PropEdgeLabel(), result );

		  // collect all nodes participating in aromatic edges
		std::set< size_t > aromaticNodes;
		for (EdgeSet::const_iterator e = aromaticEdges.begin();
				e != aromaticEdges.end(); ++e)
		{
			aromaticNodes.insert(e->first);
			aromaticNodes.insert(e->second);
		}

		  // relabel aromatic nodes
		boost::graph_traits<Molecule>::vertex_iterator vi, vi_end;
		boost::tie(vi, vi_end) = vertices( result );
		for (; vi != vi_end; ++vi) {
			  // check if the node should be aromatic
			const bool shouldBeAromatic = aromaticNodes.find(resultNodeIndex[*vi]) != aromaticNodes.end();
			  // get atom data of current label to check if currently aromatic
			std::string nodeLabel = resultNodeLabel[*vi];
			const std::string atomLabel = MoleculeUtil::getAtom(nodeLabel);
			const MoleculeUtil::AtomLabelData * atomData
				= MoleculeUtil::getAtomData(atomLabel);
			assert( atomData != NULL /*unknown atom label*/);
			  // check if relabeling is needed
			if (	((atomData->isAromatic != 0 && !shouldBeAromatic)
					||	(atomData->isAromatic == 0 && shouldBeAromatic))
				&& (MoleculeUtil::getAromaticPendant(atomLabel) != NULL))
			{
				  // relabel
				  // get aromatic/non-aromatic pendant
				const std::string* newAtom = MoleculeUtil::getAromaticPendant(atomLabel);
				assert( newAtom != NULL /*obviously we tried to relabel a non-aromatic atom*/);
				  // replace atom label within node label (can be complex label)
				nodeLabel.replace( nodeLabel.find(atomLabel), newAtom->size(), *newAtom );
				  // set new node label
				resultNodeLabel[*vi] = nodeLabel;
			}
		}

		  // relabel edges

		  // get access to edge descriptors
		boost::graph_traits<Molecule>::edge_iterator ei, e_end;
		boost::tie(ei, e_end) = boost::edges( result );
		EdgeSet arom2nonaromEdges;
		  // dummy label to set for all edges that are currently aromatic but
		  // have to be relabeled to non-aromatic
		const std::string UNCERTAIN_EDGE_LABEL = std::string("????");
		  // check all edges if they have to be relabeled
		for (; ei != e_end; ++ei) {
			Edge curEdge( resultNodeIndex[boost::source(*ei,result)]
								, resultNodeIndex[boost::target(*ei,result)] );
			if (curEdge.first > curEdge.second) {
				curEdge = Edge(curEdge.second, curEdge.first);
			}
			  // check if the edge should be aromatic
			const bool shouldBeAromatic = aromaticEdges.find(curEdge) != aromaticEdges.end();
			std::string bondLabel = resultEdgeLabel[*ei];
			const MoleculeUtil::BondLabelData * bondData
						= MoleculeUtil::getBondData( bondLabel );
			assert( bondData != NULL /*unknown bond label*/);
			  // check if aromatic relabeling needed
			if ( bondData->isAromatic == 0 && shouldBeAromatic ) {
				const std::string* newBond = MoleculeUtil::getAromaticPendant(bondLabel);
				assert( newBond != NULL /*obviously we tried to relabel a non-aromatic bond*/);
				resultEdgeLabel[*ei] = *newBond;
			} else
			  // check if non-aromatic relabeling needed
			if ( bondData->isAromatic != 0 && !shouldBeAromatic ) {
				arom2nonaromEdges.insert(curEdge);
				resultEdgeLabel[*ei] = UNCERTAIN_EDGE_LABEL;
			}
		}


		  // relabel formerly aromatic edges to non-aromatic

		if ( !arom2nonaromEdges.empty() ) {

			  // get labels for single and double bonds from MoleculeUtil
			const MoleculeUtil::BondDataMap & bondDataMap = MoleculeUtil::getBondData();
			std::string singleBondLabel="", doubleBondLabel="", aromBondLabel="";
			for (MoleculeUtil::BondDataMap::const_iterator bd=bondDataMap.begin();
					bd != bondDataMap.end(); ++bd)
			{
				if (bd->second.valence == 1 && bd->second.isAromatic == 0) {
					singleBondLabel = bd->first;
				}
				if (bd->second.valence == 2 && bd->second.isAromatic == 0) {
					doubleBondLabel = bd->first;
				}
				if (bd->second.valence == 1 && bd->second.isAromatic != 0) {
					aromBondLabel = bd->first;
				}
			}
			assert(singleBondLabel != "");
			assert(doubleBondLabel != "");
			assert(aromBondLabel != "");

			std::vector< AdjacencyData* > nodes2relabel(boost::num_vertices(result), NULL);
			std::vector< AdjacencyData* >::const_iterator nIt;

			size_t nodes2relabelCount = 0;
			for( EdgeSet::const_iterator e = arom2nonaromEdges.begin();
					e != arom2nonaromEdges.end(); ++e )
			{
				AdjacencyData *curNode = NULL;
				curNode = nodes2relabel[ e->first ];
				  // check if adjacent from node was already handled, if not do
				if (curNode == NULL) {
					boost::graph_traits<Molecule>::vertex_descriptor node
						= boost::vertex( e->first, result );
					int remainingValence = (int)MoleculeUtil::getAtomData( resultNodeLabel[node] )->valence;
					remainingValence += MoleculeUtil::getCharge( resultNodeLabel[node] );

					size_t adjAromaticEdges = 0;
					size_t openEdges = 0;
					  // iterate over all adjacent edges and reduce remaining valence
					boost::graph_traits<Molecule>::out_edge_iterator oei, oe_end;
					boost::tie(oei, oe_end) = boost::out_edges( node, result );
					for (; oei != oe_end; ++oei) {
						  // update information based on already labeled edges
						if (resultEdgeLabel[*oei] != UNCERTAIN_EDGE_LABEL) {
							remainingValence -= (int)MoleculeUtil::getBondData( resultEdgeLabel[*oei] )->valence;
							adjAromaticEdges += (MoleculeUtil::getBondData( resultEdgeLabel[*oei] )->isAromatic!=0?1:0);
						} else {
							  // count edges to be relabeled
							++openEdges;
						}
					}
					assert(adjAromaticEdges == 0 || adjAromaticEdges >= 2  /*either involved in at least one aromatic ring or not at all*/);
					  // insert new node information
					nodes2relabel[e->first] = new AdjacencyData( e->first, remainingValence, openEdges, adjAromaticEdges );
					++nodes2relabelCount;
				}
				  // check if adjacent target node was already handled, if not do
				curNode = nodes2relabel[ e->second ];
				if (curNode == NULL) {
					boost::graph_traits<Molecule>::vertex_descriptor node
						= boost::vertex( e->second, result );
					int remainingValence = (int)MoleculeUtil::getAtomData( resultNodeLabel[node] )->valence;
					remainingValence += MoleculeUtil::getCharge( resultNodeLabel[node] );

					size_t adjAromaticEdges = 0;
					size_t openEdges = 0;
					  // iterate over all adjacent edges and reduce remaining valence
					boost::graph_traits<Molecule>::out_edge_iterator oei, oe_end;
					boost::tie(oei, oe_end) = boost::out_edges( node, result );
					for (; oei != oe_end; ++oei) {
						  // update information based on already labeled edges
						if (resultEdgeLabel[*oei] != UNCERTAIN_EDGE_LABEL) {
							remainingValence -= (int)MoleculeUtil::getBondData( resultEdgeLabel[*oei] )->valence;
							adjAromaticEdges += (MoleculeUtil::getBondData( resultEdgeLabel[*oei] )->isAromatic!=0?1:0);
						} else {
							  // count edges to be relabeled
							++openEdges;
						}
					}
					assert(adjAromaticEdges == 0 || adjAromaticEdges >= 2  /*either involved in at least one aromatic ring or not at all*/);
					  // insert new node information
					nodes2relabel[e->second] = new AdjacencyData( e->second, remainingValence, openEdges, adjAromaticEdges );
					++nodes2relabelCount;
				}
			}
			  // get sorted access to the adjacent nodes
			std::vector< AdjacencyData* > nodes2relabelSorted(nodes2relabelCount,NULL);
			size_t i=0;
			for (nIt = nodes2relabel.begin(); nIt != nodes2relabel.end(); ++nIt) {
				if( *nIt != NULL) {
					  // store pointer to according AdjacencyData object
					nodes2relabelSorted[i] = *nIt;
					++i;
				}
			}
			assert( i == nodes2relabelCount );

//			std::cerr <<"\n DEBUG after init :\n";
//			for (size_t x=0; x < nodes2relabelSorted.size(); ++x) {
//				std::cerr <<"    " <<x <<" " <<nodes2relabelSorted.at(x);
//				if (nodes2relabelSorted.at(x) != NULL) {
//					std::cerr <<" "<< "id " <<nodes2relabelSorted.at(x)->nodeID <<" rem " <<nodes2relabelSorted.at(x)->remValence <<" oe " <<nodes2relabelSorted.at(x)->openEdges <<" ae " <<nodes2relabelSorted.at(x)->aromaticEdges;
//				}
//				std::cerr <<"\n";
//			}

			std::sort( nodes2relabelSorted.begin(), nodes2relabelSorted.end(), AdjacencyComp() );


//			std::cerr <<"\n DEBUG after first sort :\n";
//			for (size_t x=0; x < nodes2relabelSorted.size(); ++x) {
//				std::cerr <<"    " <<x <<" " <<nodes2relabelSorted.at(x);
//				if (nodes2relabelSorted.at(x) != NULL) {
//					std::cerr <<" "<< "id " <<nodes2relabelSorted.at(x)->nodeID <<" rem " <<nodes2relabelSorted.at(x)->remValence <<" oe " <<nodes2relabelSorted.at(x)->openEdges <<" ae " <<nodes2relabelSorted.at(x)->aromaticEdges;
//				}
//				std::cerr <<"\n";
//			}


			try {

			// iterate relabeling of front entry + sort till all done or
			// no deterministic relabeling possible
			while ( !nodes2relabelSorted.empty() && (*nodes2relabelSorted.begin())->openEdges > 0 )
			{

				AdjacencyData & curNode = *(*nodes2relabelSorted.begin());

//std::cerr <<" DEBUG curNode : id = "<<curNode.nodeID <<" val = " <<curNode.valence <<" edges = " <<curNode.edges <<" arom = " <<curNode.aromaticEdges<<"\n";
				// iterate over all uncertain edges and relabel
				boost::graph_traits<Molecule>::out_edge_iterator oei, oe_end;
				boost::tie(oei, oe_end) = boost::out_edges( boost::vertex( curNode.nodeID, result ), result );

				  // check if only single bonds left
				if ( curNode.remValence == curNode.openEdges ) {
					// find all edges and label with valence 1
					for (; oei != oe_end; ++oei) {
						if (resultEdgeLabel[*oei] == UNCERTAIN_EDGE_LABEL) {
							resultEdgeLabel[*oei] = singleBondLabel;
							// update adjacent node information
							size_t targetNode = resultNodeIndex[boost::target( *oei, result )];
							assert(nodes2relabel[ targetNode ] != NULL);
							// check if valence of target node compatible
							if ( checkValence && (nodes2relabel[ targetNode ]->openEdges == 0 || nodes2relabel[ targetNode ]->remValence == 0) ) {
								throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : tried to add single bond but target has no valence left");
							}
							nodes2relabel[ targetNode ]->openEdges--;
							nodes2relabel[ targetNode ]->remValence -= 1;
							// check if remaining valence of target node sufficient
							if ( checkValence && (nodes2relabel[ targetNode ]->openEdges > nodes2relabel[ targetNode ]->remValence) )
							{
								throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : after single bond insert not enough valence left within target node");
							}
						}
					}
				} else
				if ( curNode.openEdges == 1 ) {
					// find edge and label according to valence
					for (; oei != oe_end; ++oei) {
						if (resultEdgeLabel[*oei] == UNCERTAIN_EDGE_LABEL) {
							size_t edgeValence = 0;
							// relabel according to valence
							switch (curNode.remValence) {
							case 1 : resultEdgeLabel[*oei] = singleBondLabel;
								edgeValence = 1;
								break;
							case 2 :
								  // check if node is already involved in aromatic ring
								if (curNode.aromaticEdges >= 2) {
									  // one electron will be used by the aromatic ring
									  // thus only one left for a single bond
									resultEdgeLabel[*oei] = singleBondLabel;
									edgeValence = 1;
								} else {
									  // two electrons available for the bond
									resultEdgeLabel[*oei] = doubleBondLabel;
									edgeValence = 2;
								}
								break;
							case 3 :
								  // check if involved in aromatic ring -> error if not
								if (curNode.aromaticEdges < 2) {
									throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : remaining atom valence left for a single bond relabeling is too high (>2), no canonical labeling possible");
								}
								  // two electrons available for the bond
								resultEdgeLabel[*oei] = doubleBondLabel;
								edgeValence = 2;
								break;
							case 0 :
								throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : remaining atom valence for bond relabeling is zero, no relabeling possible");
							default:
								throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : remaining atom valence left for a single bond relabeling is too high (>2), no canonical labeling possible");
							}
							// update adjacent node information
							size_t targetNode = resultNodeIndex[boost::target( *oei, result )];
							assert(nodes2relabel[targetNode] != NULL);
							// check if valence of target node compatible
							if ( checkValence && (nodes2relabel[ targetNode ]->openEdges == 0 || nodes2relabel[ targetNode ]->remValence < edgeValence) ) {
								throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : tried to add bond but target has no valence left");
							}
							nodes2relabel[ targetNode ]->openEdges--;
							nodes2relabel[ targetNode ]->remValence -= edgeValence;
							// check if remaining valence of target node sufficient
							if ( checkValence && (nodes2relabel[ targetNode ]->openEdges > nodes2relabel[ targetNode ]->remValence) )
							{
//for (EdgeSet::const_iterator ae=aromaticEdges.begin(); ae!=aromaticEdges.end(); ++ae) {
//	std::cerr <<" " <<ae->first <<"-" <<ae->second;
//}
//for (std::vector< AdjacencyData* >::const_iterator nIt = nodes2relabelSorted.begin(); nIt != nodes2relabelSorted.end(); ++nIt, ++i) {
//	std::cerr <<" list : id = "<<(*nIt)->nodeID <<" val = " <<(*nIt)->valence <<" edges = " <<(*nIt)->edges <<" arom = " <<(*nIt)->aromaticEdges;
//}
//std::cerr	<<"\n\n"
//		<<"current relabeled graph:\n"
//		<<Molecule_Graph(result)
//		<<"\n"
//		<<"last edge = ("<<curNode.nodeID <<","<<targetNode<<")\n";

								throw std::runtime_error("ggl::chem::AromaticityPerception::relabelMolecule : after bond insert not enough valence left within target node");
							}
						}
					}
				} else {

					// NO DETERMINISTIC NON-AROMATIC LABELING POSSIBLE !!!
					// due to ordering: all following nodes have the same problem
					// or are already correctly relabeled

//					std::cout <<"\n ERROR : deterministic labeling not possible !\n cur node = "
//							<<curNode.nodeID <<" open = " <<curNode.openEdges <<" val = " <<curNode.remValence <<" arom = " <<curNode.aromaticEdges <<"\n"
//							<<Molecule_Graph(result)
//							<<std::endl;

					break;  // -> break and apply special handling

				}

				  // update curNode data
				curNode.remValence = 0;
				curNode.openEdges = 0;

				--nodes2relabelCount;

				  // update sorting due to changes in adjacent nodes' data
				std::swap( *nodes2relabelSorted.begin(), *(nodes2relabelSorted.begin()+nodes2relabelCount) );
				std::sort( nodes2relabelSorted.begin(), nodes2relabelSorted.begin()+nodes2relabelCount, AdjacencyComp() );

			}

			  // check if there are any edges left to be labeled non-aromatic
			  // but no deterministic labeling possible
			  // --> make these bonds/rings aromatic
			if (!nodes2relabelSorted.empty() && (*nodes2relabelSorted.begin())->openEdges > 0) {
				  // container that will hold all remaining edges to rename
				EdgeSet remainingEdges;
				boost::tie(ei, e_end) = boost::edges( result );
				  // check all edges if they have to be relabeled
				  // TODO slow : replace by iteration of ring bonds or dedicated container of non-aromatic bonds
				for (; ei != e_end; ++ei) {
					  // add uncertain edge to container
					if (resultEdgeLabel[*ei] == UNCERTAIN_EDGE_LABEL) {
						Edge curEdge( resultNodeIndex[boost::source(*ei,result)]
											, resultNodeIndex[boost::target(*ei,result)] );
						if (curEdge.first > curEdge.second) {
							curEdge = Edge(curEdge.second, curEdge.first);
						}
						remainingEdges.insert(curEdge);
					}
				}
				  // iterate until all edges have been assigned
				while( !remainingEdges.empty() ) {
					  // find RingDescriptor that covers most remaining edges
					RingDescriptor *nextRing = NULL;
					size_t maxEdgeOverlap = 0;
					for (size_t i=0; i<allRings.size(); ++i) {
						  // check if prediction was non-aromatic
						if (allRings.at(i)->predState == RingDescriptor::NonAromatic) {
							size_t curOverlap = 0;
							EdgeSet::const_iterator remainEdge = remainingEdges.begin(), remainEnd = remainingEdges.end();
							EdgeSet::const_iterator curEdge = allRings.at(i)->edges.begin(), curEnd = allRings.at(i)->edges.end();
							while(remainEdge != remainEnd && curEdge != curEnd) {
								  // check if overlap
								if (*remainEdge == *curEdge) {
									++curOverlap;
									++remainEdge;
									++curEdge;
								  // check which counter to increase
								} else if (*remainEdge < *curEdge) {
									++remainEdge;
								} else {
									++curEdge;
								}
							}
							  // check if new maximum found
							if (curOverlap > maxEdgeOverlap) {
								  // update maximum data
								maxEdgeOverlap = curOverlap;
								nextRing = allRings[i];
							}
						}
					}

					  // check if the remaining edges are part of a known ring
					if (nextRing == NULL) {

						// TODO check if now a deterministic labeling is possible ?!?!

						std::stringstream err;
						err	<<"\n RUNTIME ERROR in ggl::chem::AromaticityPerception::relabelMolecule :\n\n"
								<<" No deterministic non-aromatic relabeling possible. \n"
								<<"\n"
								<<" Remaining (formerly aromatic) edges that cannot be assigned are not part of any (predicted) ring.\n"
								<<" Edges = ";
						for (EdgeSet::const_iterator edge=remainingEdges.begin(); edge!=remainingEdges.end(); ++edge) {
							err <<" " <<edge->first <<"-" <<edge->second <<", ";
						}
	//					err <<"\n";
	//					for (std::vector< AdjacencyData* >::const_iterator nIt = nodes2relabelSorted.begin(); nIt != nodes2relabelSorted.end(); ++nIt, ++i) {
	//						err <<" list : id = "<<(*nIt)->nodeID <<" val = " <<(*nIt)->valence <<" edges = " <<(*nIt)->edges <<" arom = " <<(*nIt)->aromaticEdges<<"\n";
	//					}
						err	<<"\n"
								<<" Current relabeled graph:\n"
								<<Molecule_Graph(result)
								<<"\n";

						throw std::runtime_error(err.str());

					}

					  // make this ring aromatic as well
					nextRing->predState = RingDescriptor::Aromatic;

					  // relabel all nodes if not already aromatic
					for (RingList::const_iterator curID = nextRing->ring.begin(); curID != nextRing->ring.end(); ++curID) {
						  // get access to according node
						Molecule::vertex_descriptor curNode = boost::vertex( *curID, result );
						  // get atom data of current label to check if currently aromatic
						std::string nodeLabel = resultNodeLabel[curNode];
						const std::string atomLabel = MoleculeUtil::getAtom(nodeLabel);
						const MoleculeUtil::AtomLabelData * atomData
							= MoleculeUtil::getAtomData(atomLabel);
						assert( atomData != NULL /*unknown atom label*/);
						  // check if relabeling is needed
						if (	(atomData->isAromatic == 0)
							&& (MoleculeUtil::getAromaticPendant(atomLabel) != NULL))
						{
							  // relabel
							  // get aromatic/non-aromatic pendant
							const std::string* newAtom = MoleculeUtil::getAromaticPendant(atomLabel);
							assert( newAtom != NULL /*obviously we tried to relabel a non-aromatic atom*/);
							  // replace atom label within node label (can be complex label)
							nodeLabel.replace( nodeLabel.find(atomLabel), newAtom->size(), *newAtom );
							  // set new node label
							resultNodeLabel[curNode] = nodeLabel;
						}
					}

					  // relabel ring edges to aromatic labels
					for (EdgeSet::const_iterator curEdge = nextRing->edges.begin(); curEdge != nextRing->edges.end(); ++curEdge) {
						  // get access to according edge
						Molecule::edge_descriptor curResEdge = boost::edge( boost::vertex(curEdge->first,result)
																		, boost::vertex(curEdge->second,result), result ).first;
						  // check if edge label still unknown
						if (resultEdgeLabel[curResEdge] == UNCERTAIN_EDGE_LABEL) {
							  // make this an aromatic edge
							resultEdgeLabel[curResEdge] = aromBondLabel;
							  // remove from remaining edges
							remainingEdges.erase(*curEdge);
							  // skip remaining handling for this edge
							continue;
						}
						  // check if the edge should be aromatic
						std::string bondLabel = resultEdgeLabel[curResEdge];
						const MoleculeUtil::BondLabelData * bondData
									= MoleculeUtil::getBondData( bondLabel );
						assert( bondData != NULL /*unknown bond label*/);
						  // check if aromatic relabeling needed
						if ( bondData->isAromatic == 0 ) {
							resultEdgeLabel[curResEdge] = aromBondLabel;
						}
					}
				} // while remainingEdges not empty
			}

			} catch (std::runtime_error & e) {
				 // remove local data
				for (nIt = nodes2relabel.begin(); nIt != nodes2relabel.end(); ++nIt) {
					if( *nIt != NULL) {
						delete (*nIt);
					}
				}
				  // remove temporary data
				clearData();
				  // forward error to next level
				throw std::runtime_error(e.what());
			}

			 // remove local data
			for (nIt = nodes2relabel.begin(); nIt != nodes2relabel.end(); ++nIt) {
				if( *nIt != NULL) {
					delete (*nIt);
				}
			}

		} // IF relabeling aromatic -> non-aromatic necessary

	}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


	AromaticityPerception::
	RingDescriptor::
	RingDescriptor()
	 : edges()
		, ring()
		, predState(Unknown)
		, predCertainty(0)
	{
	}


////////////////////////////////////////////////////////////////////////////////


	AromaticityPerception::
	RingDescriptor::
	RingDescriptor(	const RingList & ringList )
	 : edges()
		, ring(ringList)
		, predState(Unknown)
		, predCertainty(0)
	{
		  // generate edge list description
		RingList::const_iterator cur=ringList.begin(), last=ringList.begin();
		for (cur++; cur!=ringList.end(); ++cur,++last) {
			if (*cur<*last) {
				edges.insert( Edge(*cur,*last) );
			} else {
				edges.insert( Edge(*last,*cur) );
			}
		}
	}


////////////////////////////////////////////////////////////////////////////////

}} // namespaces
