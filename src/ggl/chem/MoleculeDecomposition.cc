
#include <sstream>
#include <cstdlib>
#include <iterator>

#include "ggl/chem/MoleculeDecomposition.hh"
#include "ggl/chem/MoleculeComponent_GMLparser.hh"

#include "sgm/Match_Reporter.hh"

namespace ggl {
 namespace chem {

//////////////////////////////////////////////////////////////////////////////

	 MoleculeDecomposition::ComponentContainer
		 MoleculeDecomposition::groups = ComponentContainer();

	 MoleculeDecomposition::ComponentContainer
		 MoleculeDecomposition::interactions = ComponentContainer();

	 MoleculeDecomposition::ComponentContainer
		 MoleculeDecomposition::smallMolecules = ComponentContainer();

	 double
		 MoleculeDecomposition::correctionHydroCarbon = 3.68;
	 double
		 MoleculeDecomposition::correctionPerHeteroaromaticRing = -1.95;
	 double
		 MoleculeDecomposition::correctionPerThreeMemberedRing = 14.4;
	 double
		 MoleculeDecomposition::correctionPerAmide = -14.3;
	 double
		 MoleculeDecomposition::correctionPerThioester = -11.3;
	 double
		 MoleculeDecomposition::correctionVicinalChlorine = 1.92;

	 double
		 MoleculeDecomposition::contributionMiddlePhosphate = -208;

//////////////////////////////////////////////////////////////////////////////

MoleculeDecomposition
::RingDescriptor
::RingDescriptor( const sgm::RingReporter::RingList & ringList )
 : edges()
{
	assert(ringList.size()>1);

	  // add ring to edge set
	RingList::const_iterator node=ringList.begin(), lastNode = ringList.begin();
	for (node++; node!=ringList.end(); ++node, ++lastNode) {
		if (*node < *lastNode) {
			edges.insert(Edge(*node,*lastNode));
		} else {
			edges.insert(Edge(*lastNode,*node));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

bool
MoleculeDecomposition
::RingDescriptor
::operator<(const RingDescriptor& r2) const
{
	  // decide order on ring length (edge set as tie breaker)
	return edges.size() < r2.edges.size()
			|| (edges.size() == r2.edges.size() && edges < r2.edges);
}


//////////////////////////////////////////////////////////////////////////////

bool
MoleculeDecomposition
::RingDescriptor
::isContained( const RingDescriptor& largerRing ) const
{
	  // check if ring is small enough to be contained
	if ( edges.size() >= largerRing.edges.size() ) {
		return false;
	}
	  // check for containment with only one edge not common
	RingDescriptor::EdgeSet diffRemainder;
	std::insert_iterator< RingDescriptor::EdgeSet > insert_it (diffRemainder,diffRemainder.begin());
	std::set_difference( edges.begin(), edges.end()
					, largerRing.edges.begin(), largerRing.edges.end()
					, insert_it );
	  // if only one edge is left : ring is contained
	return diffRemainder.size() == 1;
}


//////////////////////////////////////////////////////////////////////////////

const MoleculeDecomposition::ComponentContainer &
MoleculeDecomposition
:: getGroups( void ) throw (std::logic_error)
{
	  // check if groups have to be filled
	if (groups.empty()) {

		// include the GML definitions of the MoleculeComponents to be matched
		const std::string groupsGML[] = {
			#include "ggl/chem/MoleculeDecomposition_groups.icc"
			///////////////////////////////////////////////////////////////////
			"" // THE EMPTY STRING MARKS THE END OF THE COMPONENT LIST
			///////////////////////////////////////////////////////////////////
			};

		for (size_t i=0; !(groupsGML[i].empty()); ++i ) {

			try {
				  // parse current component
				std::pair<MoleculeComponent, int> parsed
					= MoleculeComponent_GMLparser::parseGML(groupsGML[i]);
				  // ensure that parsing was successful
				if (parsed.second >= 0 ){
					throw std::logic_error("MoleculeComponent_GMLparser : "
						"parsing of group to match failed ! '"+groupsGML[i]+"'");
				}
				  // ensure priority and energy is given
				if ( parsed.first.priority != parsed.first.priority ) {
					throw std::invalid_argument("MoleculeComponent_GMLparser : "
							"no priority given");
				}
				if ( parsed.first.freeEnergy != parsed.first.freeEnergy ) {
					throw std::invalid_argument("MoleculeComponent_GMLparser : "
							"no energy given");
				}

				  // store parsed component
				groups.insert(parsed.first);

			} catch  (std::exception& ex) {
				throw std::logic_error("MoleculeDecomposition : initialization"
					" of groups to match failed ! Should never occur!"
					"\n Reported Error :\n"+std::string(ex.what()));
			}
		}
	}
	return groups;
}


//////////////////////////////////////////////////////////////////////////////


const MoleculeDecomposition::ComponentContainer &
MoleculeDecomposition
:: getInteractions( void ) throw (std::logic_error)
{
	  // check if interactions have to be filled
	if (interactions.empty()) {

		// include the GML definitions of the MoleculeComponents to be matched
		const std::string interactionsGML[] = {
#include "ggl/chem/MoleculeDecomposition_interactions.icc"
			///////////////////////////////////////////////////////////////////
			"" // THE EMPTY STRING MARKS THE END OF THE COMPONENT LIST
			///////////////////////////////////////////////////////////////////
			};

		for (size_t i=0; !(interactionsGML[i].empty()); ++i ) {

			  // parse current component
			std::pair<MoleculeComponent, int> parsed
				= MoleculeComponent_GMLparser::parseGML(interactionsGML[i]);

			  // ensure that parsing was successful
			if (parsed.second >= 0 ){
				throw std::logic_error("MoleculeDecomposition : initialization"
					" of interactions to match failed ! Should never occur!");
			}
			  // store parsed component
			interactions.insert(parsed.first);
		}
	}
	return interactions;
}


//////////////////////////////////////////////////////////////////////////////


const MoleculeDecomposition::ComponentContainer &
MoleculeDecomposition
:: getSmallMolecules( void ) throw (std::logic_error)
{
	  // check if smallMolecules have to be filled
	if (smallMolecules.empty()) {

		// include the GML definitions of the MoleculeComponents to be matched
		const std::string interactionsGML[] = {
#include "ggl/chem/MoleculeDecomposition_singleAtoms.icc"
#include "ggl/chem/MoleculeDecomposition_smallMolecules.icc"
			///////////////////////////////////////////////////////////////////
			"" // THE EMPTY STRING MARKS THE END OF THE COMPONENT LIST
			///////////////////////////////////////////////////////////////////
			};

		for (size_t i=0; !(interactionsGML[i].empty()); ++i ) {

			  // parse current component
			std::pair<MoleculeComponent, int> parsed
				= MoleculeComponent_GMLparser::parseGML(interactionsGML[i]);

			  // ensure that parsing was successful
			if (parsed.second >= 0 ){
				throw std::logic_error("MoleculeDecomposition : initialization"
					" of smallMolecules to match failed ! Should never occur!");
			}
			  // store parsed component
			smallMolecules.insert(parsed.first);
		}
	}
	return smallMolecules;
}


//////////////////////////////////////////////////////////////////////////////


size_t
MoleculeDecomposition
:: getAmideNumber( sgm::SubGraphMatching & matcher
					, sgm::Graph_Interface & molGraph )
{
	std::set< size_t > amideNitrogen;

	  // create amide pattern for matching
	Molecule amide;
	  // access data structures
	boost::property_map< Molecule, PropNodeLabel >::type
			nodeLabel(boost::get( PropNodeLabel(), amide ));
	boost::property_map< Molecule, PropEdgeLabel>::type
			edgeLabel(boost::get( PropEdgeLabel(), amide ));
	Molecule::vertex_descriptor atom;
	Molecule::vertex_descriptor atomNitrogen;
	Molecule::vertex_descriptor atomCarbon;
	Molecule::edge_descriptor edge;
	Molecule::edge_descriptor edgeNitrogen1;
	Molecule::edge_descriptor edgeNitrogen2;
	Molecule::edge_descriptor edgeNitrogenCarbon;
	Molecule::edge_descriptor edgeCarbon1;
	  // define non-aromatic pattern
	atom = boost::add_vertex( amide ); nodeLabel[atom] = "N"; // 0
		atomNitrogen = atom;
	atom = boost::add_vertex( amide ); nodeLabel[atom] = "C"; // 1
		atomCarbon = atom;
	atom = boost::add_vertex( amide ); nodeLabel[atom] = MoleculeUtil::AtomLabelWildcard; // 2
	atom = boost::add_vertex( amide ); nodeLabel[atom] = MoleculeUtil::AtomLabelWildcard; // 3
	edge = boost::add_edge( boost::vertex(0,amide), boost::vertex(1,amide), amide ).first; edgeLabel[edge] = "-";
		edgeNitrogenCarbon = edge;
	edge = boost::add_edge( boost::vertex(0,amide), boost::vertex(2,amide), amide ).first; edgeLabel[edge] = "-";
		edgeNitrogen1 = edge;
	edge = boost::add_edge( boost::vertex(0,amide), boost::vertex(3,amide), amide ).first; edgeLabel[edge] = "-";
		edgeNitrogen2 = edge;
	atom = boost::add_vertex( amide ); nodeLabel[atom] = "O"; // 4
	edge = boost::add_edge( boost::vertex(1,amide), boost::vertex(4,amide), amide ).first; edgeLabel[edge] = "=";
	atom = boost::add_vertex( amide ); nodeLabel[atom] = MoleculeUtil::AtomLabelWildcard; // 5
	edge = boost::add_edge( boost::vertex(1,amide), boost::vertex(5,amide), amide ).first; edgeLabel[edge] = "-";
		edgeCarbon1 = edge;

	  // create pattern graph object to be matched
	sgm::Graph_boost< Molecule
					, ggl::PropNodeLabel
					, ggl::PropEdgeLabel
					, ggl::PropNodeIndex> amideGraph(amide);
	sgm::Pattern pattern( amideGraph, MoleculeUtil::AtomLabelWildcard );

	// match reporter that stores nitrogen atom index of each match
	class MR_StoreIndex : public sgm::Match_Reporter {
	protected:
		std::set< size_t > & indices;
		const size_t patternNodeIndex;
	public:

		MR_StoreIndex( std::set< size_t > & indices
						, size_t patternNodeIndex )
		 : indices(indices)
			, patternNodeIndex(patternNodeIndex)
		{}

		virtual ~MR_StoreIndex()
		{}

		virtual void
		reportHit ( const sgm::Pattern_Interface & pattern,
					const sgm::Graph_Interface & target,
					const sgm::Match & match )
		{
			  // store matched index
			indices.insert(match.at(patternNodeIndex));
		}

	} mrStoreIndex( amideNitrogen, 0 );

	  // find all occurrences of the component pattern
	matcher.findMatches( pattern, molGraph, mrStoreIndex, UINT_MAX );

	  // redefine to aromatic pattern 1
	nodeLabel[atomNitrogen] = "n";
	edgeLabel[edgeNitrogen1] = ":";
	edgeLabel[edgeNitrogen2] = ":";
	  // find all occurrences of the component pattern
	matcher.findMatches( pattern, molGraph, mrStoreIndex, UINT_MAX );

	  // redefine to aromatic pattern 2
	nodeLabel[atomCarbon] = "c";
	edgeLabel[edgeCarbon1] = ":";
	edgeLabel[edgeNitrogenCarbon] = ":";
	  // find all occurrences of the component pattern
	matcher.findMatches( pattern, molGraph, mrStoreIndex, UINT_MAX );

	  // redefine to aromatic pattern 2
	edgeLabel[edgeNitrogen2] = "-";
	  // find all occurrences of the component pattern
	matcher.findMatches( pattern, molGraph, mrStoreIndex, UINT_MAX );

	return amideNitrogen.size();
}

//////////////////////////////////////////////////////////////////////////////


size_t
MoleculeDecomposition
:: getThioesterNumber( sgm::SubGraphMatching & matcher
						, sgm::Graph_Interface & molGraph )
{
	std::set< size_t > thioesterSulfur;

	  // create thioester pattern for matching
	Molecule thioester;
	  // access data structures
	boost::property_map< Molecule, PropNodeLabel >::type
			nodeLabel(boost::get( PropNodeLabel(), thioester ));
	boost::property_map< Molecule, PropEdgeLabel>::type
			edgeLabel(boost::get( PropEdgeLabel(), thioester ));
	Molecule::vertex_descriptor atom;
	Molecule::edge_descriptor edge;
	  // define pattern
	atom = boost::add_vertex( thioester ); nodeLabel[atom] = "S"; // 0
	atom = boost::add_vertex( thioester ); nodeLabel[atom] = "C"; // 1
	atom = boost::add_vertex( thioester ); nodeLabel[atom] = MoleculeUtil::AtomLabelWildcard; // 2
	edge = boost::add_edge( boost::vertex(0,thioester), boost::vertex(1,thioester), thioester ).first; edgeLabel[edge] = "-";
	edge = boost::add_edge( boost::vertex(0,thioester), boost::vertex(2,thioester), thioester ).first; edgeLabel[edge] = "-";
	atom = boost::add_vertex( thioester ); nodeLabel[atom] = "O"; // 3
	edge = boost::add_edge( boost::vertex(1,thioester), boost::vertex(3,thioester), thioester ).first; edgeLabel[edge] = "=";
	atom = boost::add_vertex( thioester ); nodeLabel[atom] = MoleculeUtil::AtomLabelWildcard; // 4
	edge = boost::add_edge( boost::vertex(1,thioester), boost::vertex(4,thioester), thioester ).first; edgeLabel[edge] = "-";

	  // create pattern graph object to be matched
	sgm::Graph_boost< Molecule
					, ggl::PropNodeLabel
					, ggl::PropEdgeLabel
					, ggl::PropNodeIndex> thioesterGraph(thioester);
	sgm::Pattern pattern( thioesterGraph, MoleculeUtil::AtomLabelWildcard );

	// match reporter that stores nitrogen atom index of each match
	class MR_StoreIndex : public sgm::Match_Reporter {
	protected:
		std::set< size_t > & indices;
		const size_t patternNodeIndex;
	public:

		MR_StoreIndex( std::set< size_t > & indices
						, size_t patternNodeIndex )
		 : indices(indices)
			, patternNodeIndex(patternNodeIndex)
		{}

		virtual ~MR_StoreIndex()
		{}

		virtual void
		reportHit ( const sgm::Pattern_Interface & pattern,
					const sgm::Graph_Interface & target,
					const sgm::Match & match )
		{
			  // store matched index
			indices.insert(match.at(patternNodeIndex));
		}

	} mrStoreIndex( thioesterSulfur, 0 );

	  // find all occurrences of the component pattern
	matcher.findMatches( pattern, molGraph, mrStoreIndex, UINT_MAX );


	return thioesterSulfur.size();
}

//////////////////////////////////////////////////////////////////////////////

void
MoleculeDecomposition
:: correctCyclicPhosphates( std::set< std::set<size_t> > & curMatchedIDs )
{
	  // get array access to matches
	std::vector< const std::set<size_t>* > matches(curMatchedIDs.size(),NULL);
	size_t i=0;
	for ( std::set< std::set<size_t> >::const_iterator s=curMatchedIDs.begin();
			s!=curMatchedIDs.end(); ++s,++i)
	{
		matches[i] = &(*s);
	}

	  // P-neighboring : what match is neighbored to what other via phosphorus
	std::vector< int > pNeigh(matches.size(),-1);
	  // O-neighboring : what match is neighbored to what other via oxygen
	std::vector< int > oNeigh(matches.size(),-1);

	  // gather neighboring information between all matches
	std::set< size_t > diffRemainder;
	std::insert_iterator< std::set< size_t > > insert_it (diffRemainder,diffRemainder.begin());
	for ( i=0; i<matches.size(); ++i ) {
		for ( size_t j=i+1; j<matches.size(); ++j) {
			diffRemainder.clear();
			std::set_difference( matches.at(i)->begin(), matches.at(i)->end()
							, matches.at(j)->begin(), matches.at(j)->end()
							, insert_it );
			if (diffRemainder.size() == 1) {
				  // phosphorus share
				pNeigh[i] = j;
				pNeigh[j] = i;
			} else if (diffRemainder.size() == (matches.at(i)->size()-1)) {
				  // oxygen share
				oNeigh[i] = j;
				oNeigh[j] = i;
			}
		}
	}

	  // generate new sorted set of matches
	std::set< std::set<size_t> > prunedMatchedIDs;
	  // handle cycle fragments
	for ( i=0; i<matches.size(); ++i ) {
		  // check if current oxygen is a fragment end
		if ( oNeigh.at(i) == -1 ) {
			  // traverse fragment
			size_t j = i;
			do {
				  // add to container
				prunedMatchedIDs.insert(*(matches.at(j)));
				  // mark as handled
				oNeigh[j] = -2;
				  // go to mirrored match and check if to continue or not
				j = pNeigh[j];
				if (oNeigh[j] >= 0) {
					size_t t = j;
					  // go to next phosphate in line
					j = oNeigh[j];
					  // mark mirrored match as handled
					oNeigh[t] = -2;
				}
			} while (oNeigh[j] >= 0);
			  // mark end as handled
			oNeigh[j] = -2;
		}
	}

	  // handle closed phosphate cycles
	for ( i=0; i<matches.size(); ++i ) {
		  // check if part of cycle -> if so start traversal here
		if (oNeigh.at(i) >= 0) {
			  // traverse cycle
			size_t j = i;
			do {
				  // add to container
				prunedMatchedIDs.insert(*(matches.at(j)));
				  // mark as handled
				oNeigh[j] = -2;
				  // go to mirrored match and check if to continue or not
				j = pNeigh[j];
				if (oNeigh[j] >= 0) {
					size_t t = j;
					  // go to next phosphate in line
					j = oNeigh[j];
					  // mark mirrored match as handled
					oNeigh[t] = -2;
				}
			} while (oNeigh[j] >= 0);
		}
	}

	  // overwrite matches
	curMatchedIDs = prunedMatchedIDs;
}

//////////////////////////////////////////////////////////////////////////////

void
MoleculeDecomposition
:: correctEnclosedPhosphates( std::set< std::set<size_t> > & curMatchedIDs )
{
	  // get array access to matches
	std::vector< const std::set<size_t>* > matches(curMatchedIDs.size(),NULL);
	size_t i=0;
	for ( std::set< std::set<size_t> >::const_iterator s=curMatchedIDs.begin();
			s!=curMatchedIDs.end(); ++s,++i)
	{
		matches[i] = &(*s);
	}

	  // neighboring : what match is neighbored to what other
	std::vector< std::vector< int > > neigh(matches.size());
	std::vector< int > neighCount(matches.size(),0);

	  // gather neighboring information between all matches
	std::set< size_t > diffRemainder;
	std::insert_iterator< std::set< size_t > > insert_it (diffRemainder,diffRemainder.begin());
	for ( i=0; i<matches.size(); ++i ) {
		for ( size_t j=i+1; j<matches.size(); ++j) {
			diffRemainder.clear();
			std::set_difference( matches.at(i)->begin(), matches.at(i)->end()
							, matches.at(j)->begin(), matches.at(j)->end()
							, insert_it );
			if (diffRemainder.size() == (matches.at(i)->size()-1)) {
				  // oxygen share
				neigh[i].push_back(j);
				neighCount[i]++;
				neigh[j].push_back(i);
				neighCount[j]++;
			}
		}
	}

	  // generate new sorted set of matches
	std::set< std::set<size_t> > prunedMatchedIDs;
	  // prune fragments to one match at one end
	for ( i=0; i<matches.size(); ++i ) {
		  // check if current match is a fragment end
		if ( neighCount.at(i) == 1 ) {
			  // add to container
			prunedMatchedIDs.insert(*(matches.at(i)));
			  // traverse fragment and clear neigh data
			int last = -1, cur = i, next = i;
			do {
				  // find successor
				for (size_t x=0; x<neigh.at(cur).size(); ++x) {
					if (neigh.at(cur).at(x) != last) {
						next = neigh.at(cur).at(x);
						break;
					}
				}
				  // mark as handled
				neigh[cur].clear();
				neighCount[cur] = -1;
				  // update setup
				last = cur;
				cur = next;
			} while (neigh.at(cur).size() > 1);
			  // mark end as handled
			neigh[cur].clear();
			neighCount[cur] = -1;
		} else if ( neighCount.at(i) == 0 ) {
			  // single phosphate -> no chain
			  // add to container
			prunedMatchedIDs.insert(*(matches.at(i)));
		}
	}

	  // overwrite matches
	curMatchedIDs = prunedMatchedIDs;
}

//////////////////////////////////////////////////////////////////////////////

void
MoleculeDecomposition
::reportRing(	const sgm::Graph_Interface& graph
				, const sgm::RingReporter::RingList & ringList )
{
	  // TODO check if only to be considered for limited ring size --> NEEDED, see test 415/420

	  // store information on what nodes are involved in rings
	curMolRingNodes.insert(ringList.begin(), ringList.end());

	  // check if ring is a heteroaromatic or aromatic hydrocarbon ring
	bool onlyAromaticBonds = true;
	bool notOnlyCarbons = false;
	RingList::const_iterator node=ringList.begin(), last=ringList.begin();
	for(node++; onlyAromaticBonds && node != ringList.end(); ++node,++last) {
		  // check if current node is no carbon if only carbons so far
		const std::string nodeLabel = MoleculeUtil::getAtom(graph.getNodeLabel( *node ));
		notOnlyCarbons = notOnlyCarbons || (nodeLabel != "C" && nodeLabel != "c");
		  // check if bond label is aromatic
		assert( graph.getEdgesBegin(*node,*last) != graph.getEdgesEnd(*node,*last) /*otherwise no edge existent*/);
		onlyAromaticBonds = MoleculeUtil::getBondData( graph.getEdgesBegin(*node,*last)->getEdgeLabel() )->isAromatic != 0;
	}
	  // set ring fragment type
	MoleculeComponent::RingFragmentType ringFragmentType = MoleculeComponent::RF_undefined;
	if (onlyAromaticBonds) {
		if (notOnlyCarbons) {
			  // ring is heteroaromatic
			ringFragmentType = MoleculeComponent::RF_heteroaromatic;
		} else {
			  // ring is aromatic hydrocarbon
			ringFragmentType = MoleculeComponent::RF_aromaticHydrocarbon;
		}
	} else {
		  // non-aromatic
		ringFragmentType = MoleculeComponent::RF_nonaromatic;
	}

	  // store ring descriptor according to ring type
	curMolRings[ringFragmentType].push_back( RingDescriptor( ringList ) );

}


//////////////////////////////////////////////////////////////////////////////



bool
MoleculeDecomposition
:: priority_order
:: operator ()(	const MoleculeComponent & c1,
				const MoleculeComponent & c2 ) const
{
	  // order by priority and use pattern size and energy as tiebreaker
	return
			c1.priority < c2.priority
			|| (c1.priority == c2.priority && c1.compIDs.size() > c2.compIDs.size())
			|| (c1.priority == c2.priority && c1.compIDs.size() == c2.compIDs.size() && c1.freeEnergy < c2.freeEnergy)
			|| (c1.priority == c2.priority && c1.compIDs.size() == c2.compIDs.size() && c1.freeEnergy == c2.freeEnergy && c1.description < c2.description);
}

//////////////////////////////////////////////////////////////////////////////



double
MoleculeDecomposition
:: getEnergy( const Molecule & mol )
{
	  // reset data
	energy = 0.0;
	curComponent = NULL;
	curMatchedIDs.clear();
	curMolRingNodes.clear();
	curMolRings.clear();

	  // graph interface object of the molecule
	Molecule_Graph molGraph(mol);

	/////////////////  FIND ALL RINGS IN MOLECULE  ///////////////////////////

	  // exhaustive ring finding algorithm
	sgm::RP_Hanser96 ringFinder;

	  // find all rings and store ring information in own data structures
	  // TODO restrict maximal ring size to be considered
	ringFinder.findRings( molGraph, *this );

	  // sort ring string representations by size to enable correct heteroaromatic ring count
	for (std::map< MoleculeComponent::RingFragmentType, std::vector< RingDescriptor > >::iterator ringVec = curMolRings.begin();
			ringVec != curMolRings.end(); ++ringVec)
	{
		std::sort(ringVec->second.begin(), ringVec->second.end());
	}


	//////////////////  PRUNE FUSED RINGS  /////////////////////////////////////

	  // prune fused non-aromatic rings
	if (curMolRings.find(MoleculeComponent::RF_nonaromatic) != curMolRings.end()) {
		  // access to vector of rings to prune
		std::vector< RingDescriptor > & toPrune = curMolRings[MoleculeComponent::RF_nonaromatic];
		  // prune from large to small
		for (int i=toPrune.size()-1; i>=0; i--) {
			bool canBePruned = false;
			  // check for contained aromatic hydrocarbon rings
			if (!canBePruned && curMolRings.find(MoleculeComponent::RF_aromaticHydrocarbon) != curMolRings.end()) {
				const std::vector< RingDescriptor > & toCheck = curMolRings.find(MoleculeComponent::RF_aromaticHydrocarbon)->second;
				for( size_t c=0; !canBePruned && c<toCheck.size(); ++c ) {
					if (toCheck.at(c).edges.size() >= toPrune.at(i).edges.size()) {
						break;
					}
					  // check for containment of this ring within toPrune(i)
					canBePruned = toCheck.at(c).isContained( toPrune.at(i) );
				}
			}
			  // check for contained heteroaromatic rings
			if (!canBePruned && curMolRings.find(MoleculeComponent::RF_heteroaromatic) != curMolRings.end()) {
				const std::vector< RingDescriptor > & toCheck = curMolRings.find(MoleculeComponent::RF_heteroaromatic)->second;
				for( size_t c=0; !canBePruned && c<toCheck.size(); ++c ) {
					if (toCheck.at(c).edges.size() >= toPrune.at(i).edges.size()) {
						break;
					}
					  // check for containment of this ring within toPrune(i)
					canBePruned = toCheck.at(c).isContained( toPrune.at(i) );
				}
			}
			  // check for contained aromatic non-aromatic rings
			if (!canBePruned) {
				const std::vector< RingDescriptor > & toCheck = curMolRings.find(MoleculeComponent::RF_nonaromatic)->second;
				for( size_t c=0; !canBePruned && c<toCheck.size(); ++c ) {
					if (toCheck.at(c).edges.size() >= toPrune.at(i).edges.size()) {
						break;
					}
					  // check for containment of this ring within toPrune(i)
					canBePruned = toCheck.at(c).isContained( toPrune.at(i) );
				}
			}
			  // prune ring if obsolete
			if (canBePruned) {
				  // remove ring descriptor at position i
				toPrune.erase( toPrune.begin() + (size_t)i );
			}

		}
	}
	  // prune fused heteroaromatic rings
	if (curMolRings.find(MoleculeComponent::RF_heteroaromatic) != curMolRings.end()) {
		  // access to vector of rings to prune
		std::vector< RingDescriptor > & toPrune = curMolRings[MoleculeComponent::RF_heteroaromatic];
		  // prune from large to small
		for (int i=toPrune.size()-1; i>=0; i--) {
			bool canBePruned = false;
			  // check for contained aromatic hydrocarbon rings
			if (!canBePruned && curMolRings.find(MoleculeComponent::RF_aromaticHydrocarbon) != curMolRings.end()) {
				const std::vector< RingDescriptor > & toCheck = curMolRings.find(MoleculeComponent::RF_aromaticHydrocarbon)->second;
				for( size_t c=0; !canBePruned && c<toCheck.size(); ++c ) {
					if (toCheck.at(c).edges.size() >= toPrune.at(i).edges.size()) {
						break;
					}
					  // check for containment of this ring within toPrune(i)
					canBePruned = toCheck.at(c).isContained( toPrune.at(i) );
				}
			}
			  // check for contained heteroaromatic rings
			if (!canBePruned) {
				const std::vector< RingDescriptor > & toCheck = curMolRings.find(MoleculeComponent::RF_heteroaromatic)->second;
				for( size_t c=0; !canBePruned && c<toCheck.size(); ++c ) {
					if (toCheck.at(c).edges.size() >= toPrune.at(i).edges.size()) {
						break;
					}
					  // check for containment of this ring within toPrune(i)
					canBePruned = toCheck.at(c).isContained( toPrune.at(i) );
				}
			}
			  // prune ring if obsolete
			if (canBePruned) {
				  // remove ring descriptor at position i
				toPrune.erase( toPrune.begin() + (size_t)i );
			}

		}
	}
	  // prune fused aromatic hydrocarbon rings
	if (curMolRings.find(MoleculeComponent::RF_aromaticHydrocarbon) != curMolRings.end()) {
		  // access to vector of rings to prune
		std::vector< RingDescriptor > & toPrune = curMolRings[MoleculeComponent::RF_aromaticHydrocarbon];
		  // prune from large to small
		for (int i=toPrune.size()-1; i>=0; i--) {
			bool canBePruned = false;
			  // check for contained aromatic hydrocarbon rings
			if (!canBePruned) {
				const std::vector< RingDescriptor > & toCheck = curMolRings.find(MoleculeComponent::RF_aromaticHydrocarbon)->second;
				for( size_t c=0; !canBePruned && c<toCheck.size(); ++c ) {
					if (toCheck.at(c).edges.size() >= toPrune.at(i).edges.size()) {
						break;
					}
					  // check for containment of this ring within toPrune(i)
					canBePruned = toCheck.at(c).isContained( toPrune.at(i) );
				}
			}
			  // prune ring if obsolete
			if (canBePruned) {
				  // remove ring descriptor at position i
				toPrune.erase( toPrune.begin() + (size_t)i );
			}

		}
	}

	//////////////////////  MATCH ALL GROUPS  ////////////////////////////////

	  // copy the molecule since we have to change it later
	Molecule m(mol);
	  // write access to all node labels
	NodeLabelMap
		mNodeLabel(boost::get( PropNodeLabel(), m ));
	boost::property_map< Molecule, PropEdgeLabel>::type
		mEdgeLabel(boost::get( PropEdgeLabel(), m ));
	NodeIndexMap
		mNodeIndex(boost::get( PropNodeIndex(), m ));

	  // will be true if only H and C atoms are present
	bool curMolIsHydroCarbon = true;
	  // holds the information which carbon atoms have adjacent chlorine
	  // atoms and their number --> to compute vicinal chlorine correction
	std::map< size_t, size_t > carbonChlorinNum;

	////////////////////////////////////////////////////////////////////
	// RELABEL NODES TO BE CONFORM WITH PATTERN
	// -> no class
	// -> no implicit protons
	// -> charge explicit, i.e. "+1" instead of just "+"
	Molecule::vertex_iterator vIt, vItEnd;
	Molecule::vertex_descriptor newNode;
	Molecule::edge_descriptor newEdge;
	for (boost::tie(vIt,vItEnd) = boost::vertices(m); vIt != vItEnd; ++vIt) {
		  // get information from node label
		const std::string& label = mNodeLabel[*vIt];
		const std::string atomLabel = MoleculeUtil::getAtom(label);
		const int charge = MoleculeUtil::getCharge( label );
		const int classID = MoleculeUtil::getClass( label );
		const size_t protons = MoleculeUtil::getProtons( label ) - (atomLabel=="H"?1:0);
		  // add explicit adjacent protons if necessary
		for (int p = protons; p > 0; p--) {
			  // add proton node
			newNode = boost::add_vertex( m );
			mNodeLabel[newNode] = "H";
			  // make adjacent
			newEdge = boost::add_edge( *vIt, newNode, m ).first;
			mEdgeLabel[newEdge] = "-";
		}
		  // relabel node if necessary
		if (charge != 0 || protons != 0 || classID != 0) {
			mNodeLabel[*vIt] = MoleculeUtil::getComplexAtomLabel( atomLabel, 0, charge, 0, true );
		}

		  // check if this molecule is a hydrocarbon
		curMolIsHydroCarbon = curMolIsHydroCarbon && (atomLabel == "C" || atomLabel == "c" || atomLabel == "H");

		  // check if chlorine
		if (atomLabel == "Cl") {
			// check for adjacent carbons and update vicinal chlorine info
			Molecule::adjacency_iterator adj, adjEnd;
			for (boost::tie(adj,adjEnd) = boost::adjacent_vertices(*vIt,m); adj != adjEnd; ++adj) {
				  // check if neighbor already known as carbon
				if (carbonChlorinNum.find(mNodeIndex[*adj]) != carbonChlorinNum.end()) {
					  // increase counter
					carbonChlorinNum[mNodeIndex[*adj]] += 1;
				} else { if (  MoleculeUtil::getAtom( mNodeLabel[*adj] ) == "C"
							|| MoleculeUtil::getAtom( mNodeLabel[*adj] ) == "c")
					{
						  // insert new counter
						carbonChlorinNum[mNodeIndex[*adj]] = 1;
					}
				}
			}
		}
	}


	  // create working copy to be matched and changed
	Molecule_Graph targetGraph(m);


	/////////////////  CHECK IF IT IS A SMALL MOLECULE  ////////////////////////

	{
		  // match reporter to check if a molecule was matched or not
		sgm::MR_Counting mr_count;
		  // get contributions for conjugation corrections
		for (ComponentContainer::const_iterator c = getSmallMolecules().begin();
				c != getSmallMolecules().end(); ++c)
		{
			  // set access to currently matched component
			curComponent = &(*c);
			  // create pattern graph object to be matched
			sgm::Graph_boost<MoleculeComponent::PatternGraph
							, ggl::PropNodeLabel
							, ggl::PropEdgeLabel
							, ggl::PropNodeIndex> patternGraph(curComponent->pattern);
			sgm::Pattern pattern(	patternGraph
									, curComponent->constraints
									, MoleculeUtil::AtomLabelWildcard);

			  // reset counter
			mr_count.resetHits();
			  // check for a match of the molecule pattern
			fullMatcher->findMatches( pattern, targetGraph, mr_count, 1 );

			  // if the molecule was matched --> handle match
			if ( mr_count.getHits() > 0 ) {
				  // set energy for the molecule match
				energy = curComponent->freeEnergy;

				  // report match
				if (reporter != NULL) {
					reporter->reportSmallMolecule( *curComponent );
					reporter->reportMatchComplete( true );
				}
			////////////////////////////////////////////////////////////////////
			/// return energy value of single molecule match
			////////////////////////////////////////////////////////////////////
				return energy;
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			}

		}
	}

	/////////////////  OTHERWISE DECOMPOSE THE MOLECULE  ///////////////////////

	// TODO INTRODUCE PRECHECKS TO AVOID UNNECESSARY MAPPINGS

	//////////////////////  MATCH ALL INTERACTIONS  ////////////////////////////

	  // get number of amides based on undecomposed (relabeled) molecule graph
	const size_t numOfAmide = getAmideNumber( *matcher, targetGraph );
	  // get number of amides based on undecomposed (relabeled) molecule graph
	const size_t numOfThioester = getThioesterNumber( *matcher, targetGraph );

	  // correction term for heteroaromatic rings
	if (curMolRings.find(MoleculeComponent::RF_heteroaromatic) != curMolRings.end()) {
		energy += (double)curMolRings.find(MoleculeComponent::RF_heteroaromatic)->second.size() * correctionPerHeteroaromaticRing;
			  // report if needed
			if (curMolRings.find(MoleculeComponent::RF_heteroaromatic)->second.size() > 0 && reporter != NULL) {
				reporter->reportInteraction( std::string("hetero aromatic rings")
											, curMolRings.find(MoleculeComponent::RF_heteroaromatic)->second.size() );
			}
	}

	  // correction for three-membered rings that is independent of the atom
	  // labels within the ring
	{
		size_t curMolNumThreeMemberedRings = 0;
		  // check for threemembered rings among all types
		typedef MoleculeComponent::RingFragmentType RFT;
		for (size_t type = 0; type < (size_t)MoleculeComponent::RF_undefined; ++type) {
			if (curMolRings.find((RFT)type) != curMolRings.end()) {
				const std::vector< RingDescriptor > & rings = curMolRings.find((RFT)type)->second;
				for(std::vector< RingDescriptor >::const_iterator ring = rings.begin();
						ring != rings.end() && ring->edges.size()<=3; ++ring)
				{
					if (ring->edges.size() == 3) {
						curMolNumThreeMemberedRings++;
					}
				}
			}
		}

		energy += (double)curMolNumThreeMemberedRings * correctionPerThreeMemberedRing;
			  // report if needed
			if (curMolNumThreeMemberedRings > 0 && reporter != NULL) {
				reporter->reportInteraction( std::string("three membered rings")
											, curMolNumThreeMemberedRings );
			}
	}

	  // check if current molecule is a hydrocarbon
	if (curMolIsHydroCarbon) {
		  // correction for hydrocarbons
		energy += correctionHydroCarbon;
		  // report if needed
		if (reporter != NULL) {
			reporter->reportInteraction( std::string("molecule is hydrocarbon")
										, 1 );
		}
	}

	  // correction for amides
	energy += (double)numOfAmide * correctionPerAmide;
		  // report if needed
		if (numOfAmide > 0 && reporter != NULL) {
			reporter->reportInteraction( std::string("amide nitrogen")
							, numOfAmide );
		}

	  // correction for thioester
	energy += (double)numOfThioester * correctionPerThioester;
		  // report if needed
		if (numOfThioester > 0 && reporter != NULL) {
			reporter->reportInteraction( std::string("thioester sulfur")
							, numOfThioester );
		}

	  // correction for vicinal chlorines
	{
		std::map< size_t, size_t >::const_iterator c1, c2;
		size_t VCLdistinct = 0;
		  // iterate over all pairs of carbons with chlorine
		for(c1 = carbonChlorinNum.begin(); c1 != carbonChlorinNum.end(); ++c1) {
			c2 = c1;
			for( c2++; c2 != carbonChlorinNum.end(); ++c2) {
				// check if adjacent in molecule
				if (boost::edge( boost::vertex(c1->first,m), boost::vertex(c2->first,m), m).second) {
					  // update correction information according to Eq.3 from the article
					VCLdistinct += std::min( c1->second, c2->second );
				}
			}
		}
		  // add final correction factor
		energy += (correctionVicinalChlorine * (double)VCLdistinct);
		  // report if needed
		if (VCLdistinct > 0 && reporter != NULL) {
			reporter->reportInteraction( std::string("vicinal chlorine distinct")
							, VCLdistinct );
		}
	}

	  // get contributions for conjugation corrections
	for (ComponentContainer::const_iterator c = getInteractions().begin();
			c != getInteractions().end(); ++c)
	{
		  // set access to currently matched component
		curComponent = &(*c);

		  // container for extended constraints for this component matching
		sgm::Pattern_Interface::ConstraintVec constraints(curComponent->constraints);

		  // create ring node constraint and make a new combined constraint container for matching
		if (!curComponent->ringFragments.empty()) {
			std::vector< bool > unconstrainedNode(boost::num_vertices(curComponent->pattern), true);
			for (MoleculeComponent::RingFragVec::const_iterator r = curComponent->ringFragments.begin();
					r != curComponent->ringFragments.end(); ++r )
			{
				for (MoleculeComponent::RingFragmentList::const_iterator n = r->fragment.begin(); n!=r->fragment.end(); ++n) {
					if (unconstrainedNode.at(*n)) {
						  // add new ring constraint for the according node
						constraints.push_back(new MC_MC_RingNode( *n, curMolRingNodes));
						  // mark that a ring node constraint was set for this node
						unconstrainedNode[*n] = false;
					}
				}
			}
		}

		  // create pattern graph object to be matched
		sgm::Graph_boost<MoleculeComponent::PatternGraph
						, ggl::PropNodeLabel
						, ggl::PropEdgeLabel
						, ggl::PropNodeIndex> patternGraph(curComponent->pattern);
		sgm::Pattern pattern(	patternGraph
								, constraints
								, MoleculeUtil::AtomLabelWildcard );

		  // find all occurrences of the component pattern
		matcher->findMatches( pattern, targetGraph, *this, UINT_MAX );

		  // if the pattern was matched --> handle matches
		if ( ! curMatchedIDs.empty() ) {
			  // increase energy for each unique match
			energy += (curComponent->freeEnergy * curMatchedIDs.size());
			  // report if needed
			if (reporter != NULL) {
				reporter->reportInteraction( curComponent->description
								, curMatchedIDs.size() );
			}
			  // reset occurrence storage
			curMatchedIDs.clear();
		}

		  // clear additional RingNode constraints
		for (size_t i=curComponent->constraints.size(); i<constraints.size();++i) {
			delete constraints[i];
		}
		constraints.clear();

	} // end for all conjugations

	////////////////////////////////////////////////////////////////////
	// DO GROUP MATCHING

	// TODO CHECK IF PROTON FILLED VERSIONS OF THE ION FORMS MAKE SENSE --> OTHERWISE DELETE


	size_t compCount = 0;
	for (ComponentContainer::const_iterator c = getGroups().begin();
			c != getGroups().end(); ++c)
	{
		  // set access to currently matched component
		curComponent = &(*c);
		  // check if the pattern requires rings and if any are available
		if ( curMolRingNodes.empty() && !curComponent->ringFragments.empty() ) {
			  // no rings present, so pattern cannot be matched .. skip it!
			continue;
		}

		  // container for extended constraints for this component matching
		sgm::Pattern_Interface::ConstraintVec constraints(curComponent->constraints);

		  // create ring node constraint and make a new combined constraint container for matching
		if (!curComponent->ringFragments.empty()) {
			std::vector< bool > unconstrainedNode(boost::num_vertices(curComponent->pattern), true);
			for (MoleculeComponent::RingFragVec::const_iterator r = curComponent->ringFragments.begin();
					r != curComponent->ringFragments.end(); ++r )
			{
				for (MoleculeComponent::RingFragmentList::const_iterator n = r->fragment.begin(); n!=r->fragment.end(); ++n) {
					if (unconstrainedNode.at(*n)) {
						  // add new ring constraint for the according node
						constraints.push_back(new MC_MC_RingNode( *n, curMolRingNodes));
						  // mark that a ring node constraint was set for this node
						unconstrainedNode[*n] = false;
					}
				}
			}
		}


		  // create pattern graph object to be matched
		sgm::Graph_boost< MoleculeComponent::PatternGraph
						, ggl::PropNodeLabel
						, ggl::PropEdgeLabel
						, ggl::PropNodeIndex> patternGraph(curComponent->pattern);
		  // create pattern graph object to be matched
		sgm::Pattern pattern(patternGraph, constraints, MoleculeUtil::AtomLabelWildcard);
		  // find all occurrences of the component pattern
		matcher->findMatches( pattern, targetGraph, *this, UINT_MAX );

		  // clear additional RingNode constraints
		for (size_t i=curComponent->constraints.size(); i<constraints.size();++i) {
			delete constraints[i];
		}
		constraints.clear();


		  // if the pattern was matched --> handle matches
		if ( ! curMatchedIDs.empty() ) {

			  // handle cyclic polyphosphates
			if (curComponent->description.find("[cyclic phosphate]") != std::string::npos) {
				  // due to symmetric pattern doubled matches that have to be pruned
				correctCyclicPhosphates( curMatchedIDs );
			}

			 // handle enclosed phosphate fragments
			if (curComponent->description.find("[primary phosphate] [enclosed]") != std::string::npos) {
				  // due to symmetric pattern the whole fragment of enclosed phosphates is matched
				  // has to be pruned to one per chain --> remaining will be covered by "middle phosphate" handling
				correctEnclosedPhosphates( curMatchedIDs );
			}

			  // increase energy for each unique match
			energy += (curComponent->freeEnergy * curMatchedIDs.size());
			  // report component matches if needed
			if (reporter != NULL) {
				reporter->reportComponent( *curComponent, curMatchedIDs.size(), compCount );
			}
			  // create set of all nodes to be relabeled to avoid double handling
			std::set<size_t> toRelabel;
			for (std::set<std::set<size_t> >::const_iterator
					s = curMatchedIDs.begin(); s != curMatchedIDs.end(); ++s )
			{
				toRelabel.insert(s->begin(), s->end());
			}
			  // relabel nodes in m for all occurrences counted
			for (std::set<size_t>::const_iterator n = toRelabel.begin();
					n != toRelabel.end(); ++n )
			{
				std::string curLabel = mNodeLabel[ boost::vertex( *n, m ) ];
				  // check if already complex atom label with class information
				if (curLabel.find(':') != std::string::npos) {
					throw std::invalid_argument("MoleculeDecomposition : "
							"molecule atom " + boost::lexical_cast<std::string>(*n) +
							" shows already a class information!"
							" Don't want to overwrite with matched component number.");
				}
				  // append component number as class information
				curLabel = curLabel + ":" + boost::lexical_cast<std::string>(compCount);
				  // store new node label
				mNodeLabel[ boost::vertex( *n, m ) ] = curLabel;
			}

			  // increased matched component counter
			compCount++;

			  // handle non-cyclic polyphosphates
			if (curComponent->description.find("[primary phosphate]") != std::string::npos) {

				size_t numOfExtensions = 0;

				  // do extension of polyphosphate starting from each matched primary phosphate
				for (std::set< std::set<size_t> >::const_iterator
						s = curMatchedIDs.begin(); s != curMatchedIDs.end(); ++s )
				{
					  // find the phosphate connecting oxygen node index (P-[O]-P)
					  // --> has to have a non-class "P" labeled neighbor
					std::set<size_t>::const_iterator compNode = s->begin();
					for (; compNode != s->end(); ++compNode) {
						  // go to node that has a "P" neighbor
						std::pair<Molecule::adjacency_iterator, Molecule::adjacency_iterator>
							neigh = boost::adjacent_vertices(boost::vertex( *compNode, m ), m );
						while(neigh.first != neigh.second) {
							if (mNodeLabel[ *(neigh.first) ] == "P") {
								  // P-neigbor found
								  // extend and relabel recursively starting at linkO
								numOfExtensions += extendPolyphosphateChain( m
																, mNodeLabel
																, mNodeIndex
																, mNodeIndex[*(neigh.first)]
																, compCount
															);
								  // stop here, should be no second unmatched P-neighbor
								break;
							}
							neigh.first++;
						}
						if (neigh.first != neigh.second) {
							  // P-neighbor found -> stop neighbor screen
							break;
						}
					}
				}

				if (numOfExtensions > 0) {
					  // update energy
					energy += numOfExtensions * contributionMiddlePhosphate;
					  // inform reporter
					if (reporter != NULL) {
						reporter->reportPolyPhosphate( numOfExtensions, compCount );
						  // increased matched component counter
						compCount++;
					}
				}
			}
			  // reset occurrence storage
			curMatchedIDs.clear();
		}


	} // end for all groups

	  // report final decomposed graph if needed
	if (reporter != NULL) {
		reporter->reportMatchGraph( targetGraph );

		  // check if all nodes were successfully matched
		bool completeMatch = true;
		for (size_t i=targetGraph.getNodeNumber(); completeMatch && i>0; i--) {
			  // check for class separator in all node labels
			completeMatch = targetGraph.getNodeLabel(i-1).find(":") != std::string::npos;
		}
		  // report match completeness
		reporter->reportMatchComplete(completeMatch);
	}

	return energy;
}


//////////////////////////////////////////////////////////////////////////////



size_t
MoleculeDecomposition
:: extendPolyphosphateChain( Molecule & m
					, MoleculeDecomposition::NodeLabelMap & mNodeLabel
					, MoleculeDecomposition::NodeIndexMap & mNodeIndex
					, const size_t idxNextP
					, const size_t compID
					)
{
	Molecule::vertex_descriptor curP = boost::vertex( idxNextP, m );

	assert( mNodeLabel[curP] == "P" /*has to be a phosphorus without label*/ );

	  // collect current phosphate nodes, i.e. first the neighbored O
	std::vector< Molecule::vertex_descriptor > neighO;
	std::pair<Molecule::adjacency_iterator, Molecule::adjacency_iterator>
		neighOfP = boost::adjacent_vertices( curP, m );
	while(neighOfP.first != neighOfP.second) {
		std::string & neighLabel = mNodeLabel[ *(neighOfP.first) ];
		if ( neighLabel.find(':') == std::string::npos ) {
			  // check if neighbor is an oxygen
			if ( MoleculeUtil::getAtom(neighLabel) == "O" ) {
				  // store neighbor
				neighO.push_back(*(neighOfP.first));
			} else {
				  // if not --> this is no phosphate --> abort
				return 0;
			}
		}
		neighOfP.first++;
	}
	  // check if complete middle phosphate matched
	if (neighO.size() != 3) {
		  // a middle phosphate has to have three unlabeled oxygens
		  // --> otherwise incomplete --> abort
		return 0;
	}

	  // relabel phoshorus
	mNodeLabel[ curP ] = MoleculeUtil::getComplexAtomLabel( "P", 0, 0, compID, true );

	  // relabel oxygens and check for adjacent phosphate
	std::vector< Molecule::vertex_descriptor > neighH;
	Molecule::vertex_descriptor successiveP = curP;
	for (size_t i=0; i<neighO.size(); ++i ) {
		std::string & neighLabel = mNodeLabel[ neighO.at(i) ];
		  // relabel neighbored oxygen
		const int charge = MoleculeUtil::getCharge(neighLabel);
		mNodeLabel[ neighO.at(i) ] = MoleculeUtil::getComplexAtomLabel( "O", 0, charge, compID, true );
		  // get neighborhood of oxygen
		std::pair<Molecule::adjacency_iterator, Molecule::adjacency_iterator>
			neighOfO = boost::adjacent_vertices( neighO.at(i), m );
		while(neighOfO.first != neighOfO.second) {
			std::string & neighLabel = mNodeLabel[ *(neighOfO.first) ];
			if (neighLabel == "H") {
				neighH.push_back( *(neighOfO.first) );
			} else if (neighLabel == "P") {
				if ( mNodeIndex[ *(neighOfO.first)] != idxNextP ) {
					successiveP = *(neighOfO.first);
				}
			}
			neighOfO.first++;
		}
	}
	  // relabel protons
	for (size_t i=0; i<neighH.size(); ++i ) {
		mNodeLabel[ neighH.at(i) ] = MoleculeUtil::getComplexAtomLabel( "H", 0, 0, compID, true );
	}

	  // check if recursion needed
	if ( successiveP == curP ) {
		  // was the last phosphate in chain
		return 1;
	} else {
		  // recursive call --> go to next phosphate in chain
		return 1 + extendPolyphosphateChain( m, mNodeLabel, mNodeIndex, mNodeIndex[successiveP], compID);
	}
}





//////////////////////////////////////////////////////////////////////////////



void
MoleculeDecomposition
:: reportHit (	const sgm::Pattern_Interface & componentPattern,
				const sgm::Graph_Interface & targetMol,
				const sgm::Match & match )
{
	//TODO DEBUG COMMENT : additional check should be redundant and thus obsolete
//	  // check if all additional matching constraints are fulfilled
//	for (MoleculeComponent::ConstraintVec::const_iterator c = curComponent->constraints.begin();
//			c != curComponent->constraints.end(); ++c )
//	{
//		  // check if current constraint is violated
//		if ( ! (*c)->isValidMatch( componentPattern, targetMol, match ) ) {
//			return;
//		}
//	}

	// TODO add constraint that node has to participate in ring to pattern constraints for faster matching

	  // check if ring constraints are fulfilled
	for (MoleculeComponent::RingFragVec::const_iterator r = curComponent->ringFragments.begin();
			r != curComponent->ringFragments.end(); ++r )
	{
		 // check current ring fragment if valid
		if ( ! isValidRingFragment(*r, match) ) {
			return;
		}
	}
	  // the set of matched component node IDs
	std::set< size_t > matchedCompIDs;
	for (MoleculeComponent::NodeSet::const_iterator n = curComponent->compIDs.begin();
			n != curComponent->compIDs.end(); ++n )
	{
		  // store current node match
		matchedCompIDs.insert( match.at(*n) );
	}
	  // store matched indices for these nodes into curMatchedIDs
	curMatchedIDs.insert( matchedCompIDs );
}


//////////////////////////////////////////////////////////////////////////////


bool
MoleculeDecomposition::
isValidRingFragment( const MoleculeComponent::RingFragment & r
					, const sgm::Match & match )
{
	  // create edge set of the ring fragment
	RingDescriptor::EdgeSet fragEdges;
	  // parallel : check if nodes of ring fragment are within rings
	  // --> this is a quick check before searching the fragment
	MoleculeComponent::RingFragmentList::const_iterator n = r.fragment.begin(), last = r.fragment.begin();
	for (n++; n != r.fragment.end(); ++n, ++last) {
		if ( curMolRingNodes.find( match.at(*n) ) == curMolRingNodes.end() ) {
			  // the current ring fragment node does not participate in
			  // any ring within the target molecule
			return false;
		}
		  // store edge for later fragment search
		if (match.at(*last) < match.at(*n)) {
			fragEdges.insert( RingDescriptor::Edge( match.at(*last),match.at(*n)) );
		} else {
			fragEdges.insert( RingDescriptor::Edge( match.at(*n),match.at(*last)) );
		}
	}

	  // get ring set description according to the fragment type
	std::map< MoleculeComponent::RingFragmentType, std::vector< RingDescriptor > >::const_iterator
		ringsInMolecule = curMolRings.find(r.type);
	 // check if there is any ring of this type
	if (ringsInMolecule == curMolRings.end()) {
		  // no ring of this type found, so abort further checks
		return false;
	}

	  // temporary variables for set difference checks
	RingDescriptor::EdgeSet diffRemainder;
	std::insert_iterator< RingDescriptor::EdgeSet > insert_it (diffRemainder,diffRemainder.begin());
	  // check for the fragment in all ring string representations
	for (size_t i=0; i<ringsInMolecule->second.size(); ++i ) {
		diffRemainder.clear();
		  // check fragment containment via set difference of edge sets
		std::set_difference( fragEdges.begin(), fragEdges.end()
						, ringsInMolecule->second.at(i).edges.begin(), ringsInMolecule->second.at(i).edges.end()
						, insert_it );
		  // if no edge is left --> all edges part of the ring --> MATCH
		if (diffRemainder.empty()) {
			return true;
		}

	}

	  // no match found --> so fragment is not part of the match
	return false;
}

//////////////////////////////////////////////////////////////////////////////


 } // namespace chem
} // namespace ggl
