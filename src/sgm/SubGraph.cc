#include "sgm/SubGraph.hh"

#include <cassert>
#include <algorithm>
#include <climits>

namespace sgm
{
////////////////////////////////////////////////////////////////////////////
	
	const SubGraph::IndexType SubGraph::NOT_MAPPED_INDEX = UINT_MAX;

////////////////////////////////////////////////////////////////////////////
		
	SubGraph
	::SubGraph(	const Graph_Interface & fullGraph_
				, const NodeList & nodeList_ )
	 :	fullGraph(fullGraph_)
	 	, nodeList(nodeList_)
	 	, full2sub( fullGraph_.getNodeNumber(), NOT_MAPPED_INDEX )
	{
		assert(nodeList.size() <= fullGraph.getNodeNumber() /* subgraph cannot be larger than source graph */);
		
		size_t srcNodeNbr = fullGraph.getNodeNumber();
		size_t localIndex = 0;
		  // calculate mapping of global to local indices
		for (NodeList::const_iterator it=nodeList.begin(); it!=nodeList.end(); ++it) {
			  // check if node indices of the subgraph are valid for the source graph
			assert(*it < srcNodeNbr /* or source node index out of bound */);
			full2sub[*it] = localIndex;
			localIndex++;
		}

	}
	
////////////////////////////////////////////////////////////////////////////
		
	SubGraph
	::SubGraph(	const Graph_Interface & fullGraph_
				, const CompLabel & compLabel
				, const size_t subGraphLabel )
	 :	fullGraph( fullGraph_ )
	 	, nodeList()
	 	, full2sub( fullGraph_.getNodeNumber(), NOT_MAPPED_INDEX )
	{
		  // for all nodes check their coloring
		for (IndexType n=0; n<compLabel.size(); ++n) {
			  // check if this node has the correct label
			if (compLabel.at((size_t)n)==subGraphLabel) {
				  // add node to set of local nodes
				nodeList.push_back(n);
				  // set up global to local mapping
				full2sub[n] = (nodeList.size()-1);
			}
		}
		
		assert(nodeList.size() > 0 /* subgraph should contain at least one node */);
	}

	

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




	SubGraph::
	EdgeDescriptor::
	EdgeDescriptor(	const Graph_Interface::OutEdge_iterator& cur_edge_,
					const Graph_Interface::OutEdge_iterator& list_end_,
					const SubGraph& parent_ )
	 : Graph_Interface::EdgeDescriptor()
		, curEdge( cur_edge_ )
		, edgeEnd( list_end_ )
		, parent(parent_)
	{
		  // go to first edge part of the subgraph
		while( curEdge != edgeEnd && parent.full2sub[curEdge->getToIndex()]==SubGraph::NOT_MAPPED_INDEX) {
			++curEdge;
		}
		  // check if no neighbors available
		if (curEdge != edgeEnd) {
			  // copy edge information
			from  = parent.full2sub[curEdge->getFromIndex()];
			to    = parent.full2sub[curEdge->getToIndex()];
		}
	}

	//////////////////////////////////////////////////////////////////////////



	SubGraph::
	EdgeDescriptor&
	SubGraph::
	EdgeDescriptor::
	operator++()
	{
		  // check if end of list already reached
		if (curEdge == edgeEnd) {
			return *this;
		}
		  // proceed to next edge part of the subgraph
		++curEdge;
		while( curEdge != edgeEnd && parent.full2sub[curEdge->getToIndex()]==SubGraph::NOT_MAPPED_INDEX) {
			++curEdge;
		}
		  // update data
		  // check if no neighbors available --> set to NO_FURTHER_EDGES
		if (curEdge != edgeEnd) {
			  // copy edge information
			to    = parent.full2sub[curEdge->getToIndex()];
		}

		return *this;
	}

	//////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


} // namespace sgm




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////  SUBGRAPHPATTERN
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

namespace sgm
{
////////////////////////////////////////////////////////////////////////////

	void
	SubGraphPattern::
	remapConstraints( const ConstraintVec & originalConstraints
						, const Match & full2sub
						, ConstraintVec & remappedConstraints )
	{
		  // clear container
		for( size_t i=0; i<remappedConstraints.size(); ++i) {
			delete(remappedConstraints[i]);
		}
		  // resize to the needed size
		remappedConstraints.resize(originalConstraints.size(), NULL);

		  // refill return container
		for( size_t i=0; i<remappedConstraints.size(); ++i) {
			remappedConstraints[i] = originalConstraints.at(i)->remap( full2sub, NOT_MAPPED_INDEX );
		}
		  // prune NULL entries
		size_t fillPos = 0;
		for( size_t i=0; i<remappedConstraints.size(); ++i) {
			  // skip empty entry
			if (remappedConstraints.at(i) == NULL) {
				continue;
			}
			  // shift non-empty entry if needed
			if (i != fillPos) {
				remappedConstraints[fillPos] = remappedConstraints[i];
			}
			++fillPos;
		}
		  // shrink container size
		remappedConstraints.resize(fillPos);
	}

////////////////////////////////////////////////////////////////////////////
} // namespace sgm
