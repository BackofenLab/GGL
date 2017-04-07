
#include "sgm/GM_vf2.hh"

 // VF2 library headers

#include <vf2/match.h>
#include <vf2/vf2_state.h>

#include <cassert>

namespace sgm {


////////////////////////////////////////////////////////////////////////////////
	// GM_vf2 class
////////////////////////////////////////////////////////////////////////////////

	size_t
	GM_vf2
	::findMatches (	const Pattern_Interface & pattern,
				const Graph_Interface & target,
				Match_Reporter & reporter,
				const size_t maxHits )
	{
		  // pre check if a graph matching is possible based on the node number
		if (pattern.getPatternGraph().getNodeNumber() != target.getNodeNumber())
			return 0;

		  // setup pattern interface
		PatternVec patterns(1, &pattern);

		  // setup reporter interface
		ReporterVec reporters(patterns.size(), &reporter);

		return VF2_MatchingHandler::findMatches
				< vf2::VF2State, NodeComparator, EdgeLabelComparator >
					(patterns, target, reporters, maxHits );

	}

////////////////////////////////////////////////////////////////////////////////

	size_t
	GM_vf2
	::findMatches (	const std::vector< const Pattern_Interface*> & inPatterns,
				const Graph_Interface & target,
				Match_Reporter & reporter,
				const size_t maxHits )
	{

		  // setup reporter interface
		PatternVec patterns(inPatterns.size(), NULL);
		for (size_t i=0; i<inPatterns.size(); ++i) {
			  // copy only if pattern and target are of equal size
			if (inPatterns.at(i)->getPatternGraph().getNodeNumber() == target.getNodeNumber())
				patterns[i] = inPatterns.at(i);
		}

		  // setup reporter interface
		ReporterVec reporters(patterns.size(), &reporter);

		return VF2_MatchingHandler::findMatches
				< vf2::VF2State, NodeComparator, EdgeLabelComparator >
					(patterns, target, reporters, maxHits );

	}

////////////////////////////////////////////////////////////////////////////////

	size_t
	GM_vf2
	::findMatches (	const std::vector< const Pattern_Interface*> & inPatterns,
				const Graph_Interface & target,
				std::vector<Match_Reporter*> & reporters,
				const size_t maxHits )
	{

		  // setup reporter interface
		PatternVec patterns(inPatterns.size(), NULL);
		for (size_t i=0; i<inPatterns.size(); ++i) {
			  // copy only if pattern and target are of equal size
			if (inPatterns.at(i)->getPatternGraph().getNodeNumber() == target.getNodeNumber())
				patterns[i] = inPatterns.at(i);
		}

		return VF2_MatchingHandler::findMatches
				< vf2::VF2State, NodeComparator, EdgeLabelComparator >
					(patterns, target, reporters, maxHits );


	}
	
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


	GM_vf2
	::NodeComparator
	::NodeComparator( const Label & pWildcard
			, const sgm::Pattern_Interface & pattern
			, const sgm::Graph_Interface & target )
	 :	LabelComparator(pWildcard)
		, pattern(pattern)
		, target(target)
	{}


////////////////////////////////////////////////////////////////////////////////

	bool
	GM_vf2
	::NodeComparator
	::compatible(void *pa, void *pb)
	{

		assert( pa != NULL /*otherwise no pattern node data available*/ );
		assert( pb != NULL /*otherwise no target node data available*/ );

		  // do cast
		NodeData* a = static_cast<NodeData*>(pa);
		NodeData* b = static_cast<NodeData*>(pb);

		  // check if node degree and labels are compatible
		bool compatible =  (a->outDegree == b->outDegree)
							&& (a->selfloops == b->selfloops)
							&& LabelComparator::compatibleLabel( a->label, b->label );
		  // check node constraints
		if (compatible && a->nodeConstraints != NULL) {
			for( NodeData::NodeConstraints::const_iterator c = a->nodeConstraints->begin();
					compatible && c != a->nodeConstraints->end(); ++c)
			{
				  // check if current constraint is fulfilled
				compatible = (*c)->isValidMatch( pattern, target, b->id );
			}
		}

		return compatible;
	}

////////////////////////////////////////////////////////////////////////////////


	GM_vf2
	::EdgeLabelComparator
	::EdgeLabelComparator( const Label& pWildcard )
	 :	LabelComparator(pWildcard)
	{}


////////////////////////////////////////////////////////////////////////////////

	bool
	GM_vf2
	::EdgeLabelComparator
	::compatible(void *pa, void *pb) {

		  // do cast
		EdgeLabel* a = static_cast<EdgeLabel*>(pa);
		EdgeLabel* b = static_cast<EdgeLabel*>(pb);

		  // ensure that label sets are of equal size
		if (a->size() != b->size()){
			return false;
		}

		  // special handling for non-multiple pattern edge
		if (a->size() == 1) {
			  // check if edge labels are compatible
			return (LabelComparator::compatibleLabel( *(a->begin()), *(b->begin()) ));
		}

		 // more complex check for multiple parallel edges
		bool allFound = true;

		EdgeLabel::const_iterator ai = a->begin();
		EdgeLabel::const_iterator aiEnd = ai;
		  // find first position that is different to *ai
		while(aiEnd != a->end() && *ai==*aiEnd) {++aiEnd;};
		  // screen through all blocks
		while( allFound && ai != aiEnd ) {
			assert( *ai != 0 );
			  // try to locate current block
			  // --> wildcard blocks are omitted, label set size check is sufficient
			if ( pWildcard == 0 || *ai != pWildcard ) {
				EdgeLabel::const_iterator bi = b->find(*ai);
				allFound = bi != b->end();
				while( allFound && ai != aiEnd ) {
					allFound = bi != b->end() && LabelComparator::compatibleLabel( *ai, *bi);
					++ai; ++bi;
				}
			}
			  // proceed to next block
			ai = aiEnd;
			if (allFound && ai != a->end()) {
				  // find first position that is different to *ai
				while(aiEnd != a->end() && *ai==*aiEnd) {++aiEnd;};
			}
		}
		  // return search result
		return allFound;
	}


////////////////////////////////////////////////////////////////////////////


} // namespace sgm

