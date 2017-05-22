
#include "MR_UniqueAtomMatch.hh"

#include <ggl/chem/Molecule.hh>
#include <ggl/chem/Molecule_Graph_noClass.hh>


//////////////////////////////////////////////////////////////////////////////

MR_UniqueAtomMatch::
MR_UniqueAtomMatch( sgm::Match_Reporter & forwardReporter_ )
 :	forwardReporter(forwardReporter_)
{
}

//////////////////////////////////////////////////////////////////////////////

MR_UniqueAtomMatch::
~MR_UniqueAtomMatch()
{
}

//////////////////////////////////////////////////////////////////////////////

void
MR_UniqueAtomMatch::
reportHit( const sgm::Pattern_Interface & pattern
			, const sgm::Graph_Interface & target
			, const sgm::Match & match )
{
	// check if target is an instance
	const ggl::chem::Molecule_Graph_V * molecules =
			dynamic_cast< const ggl::chem::Molecule_Graph_V * >( &target );

	// if not, check if wrapped
	if (molecules == NULL) {
		// check if target is an instance of an expected wrapper
		const ggl::chem::Molecule_Graph_noClass * molsNoClass =
			dynamic_cast< const ggl::chem::Molecule_Graph_noClass * >( &target );
		// if so, check wrapped graph if it can be casted
		if (molsNoClass != NULL) {
			// try cast on wrapped graph
			molecules =
				dynamic_cast< const ggl::chem::Molecule_Graph_V * >( &(molsNoClass->getWithFullAtomLabels()) );
		}
	}

	// check if atoms are uniquely mapped
	if (molecules != NULL) {

		std::set< ggl::chem::Molecule_Graph_V::LocalIndex > localIndices;
		// check if local index of each matched atom is unique
		for (size_t i=0; i<match.size(); i++) {
			// get local index = (graphID,graphNodeID)
			ggl::chem::Molecule_Graph_V::LocalIndex li = molecules->getLocalIndex( match.at(i) );

//			std::cout <<"### match "<<i <<" to "<<match.at(i)<<" = "<<molecules->getNodeLabel(match.at(i)) <<" : "<<li.first<<","<<li.second <<" : "<<(size_t)molecules->getGraphs().at(li.first)<<std::endl;
			// replace graph number with according address of the object
			li.first = (size_t)molecules->getGraphs().at(li.first);
			// check via insertion if not already present in set
			if ( ! localIndices.insert( li ).second ) {
				// stop insertion and do NOT forward match
				// since non-unique atom mapping
				return;
			}
		}
		assert( localIndices.size() == match.size() );
	}

	// forward match report since match passed all checks
	forwardReporter.reportHit( pattern, target, match );

}


//////////////////////////////////////////////////////////////////////////////
