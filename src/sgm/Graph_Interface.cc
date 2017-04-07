
#include "sgm/Graph_Interface.hh"

#include <stack>
#include <cassert>
#include <climits>
#include <set>


namespace sgm {

	
////////////////////////////////////////////////////////////////////////////

	  // Equality comparison
	  // @param toCompare the Graph_Interface to compare to
	  // @return true if both interfaces describe the same graph
	bool
	Graph_Interface
	::operator==(const Graph_Interface& toCompare ) const
	{
		  // check node numbers
		if (this->getNodeNumber() != toCompare.getNodeNumber()) {
			return false;
		}

		  // check all node labels
		for (size_t n=0; n<this->getNodeNumber(); ++n) {
			if (this->getNodeLabel(n).compare(toCompare.getNodeLabel(n))!=0) {
				return false;
			}
		}

		  // check for each node the set and order of adjacent edges
		for (size_t n=0; n<this->getNodeNumber(); ++n) {

			OutEdge_iterator thisEdge = this->getOutEdgesBegin(n);
			OutEdge_iterator thisEdgeEnd = this->getOutEdgesEnd(n);
			OutEdge_iterator compareEdge = toCompare.getOutEdgesBegin(n);
			OutEdge_iterator compareEdgeEnd = toCompare.getOutEdgesEnd(n);

			while ( thisEdge != thisEdgeEnd && compareEdge != compareEdgeEnd ) {
				  // check if edge descriptors are equal
				if (*thisEdge != *compareEdge) {
					return false;
				}
				  // go to next adjacent edge
				++thisEdge;
				++compareEdge;
			}
			  // check if no additional edges are present
			if (thisEdge != thisEdgeEnd || compareEdge != compareEdgeEnd) {
				return false;
			}
		}

		return true;
	}

	
////////////////////////////////////////////////////////////////////////////
	
	
	
	  // Inequality comparison
	  // @param toCompare the Graph_Interface to compare to
	  // @return true if both interfaces describe different graphs
	bool 
	Graph_Interface
	::operator!=(const Graph_Interface& ed ) const
	{
		return !(this->operator==(ed));
	}

////////////////////////////////////////////////////////////////////////////
	
	size_t 
	Graph_Interface
	::connectedComponents(	const Graph_Interface& g
							, Graph_Interface::CompLabel & compID ) 
	{
		const size_t nodeNumber = g.getNodeNumber();
		  // initializing the labels
		compID.resize(nodeNumber, UINT_MAX); // resizing to neccessary size
		
		  // check if anything to do
		if (nodeNumber == 0)
			return 0;
		
		  // temp data
		size_t curLabel = 0;
		size_t lastInitID = 0;
		
		// recursion for complete graph exploration
		do {
			  // color the connected component adjacent to node lastInitID
			labelAdjacentNodes( g, lastInitID, compID, curLabel);
			  // shift lastInitID to next unseen node and therefore next component
			for (lastInitID++; lastInitID < nodeNumber && compID[lastInitID] != UINT_MAX; ++lastInitID);
			  // increase connected component label
			curLabel++;
		} while (lastInitID < nodeNumber);
		
		return curLabel;
	}
	
////////////////////////////////////////////////////////////////////////////
	
	// NON-RECURSIVE VERSION
	void 
	Graph_Interface
	::labelAdjacentNodes(	const Graph_Interface& g
							, const size_t curNode
							, Graph_Interface::CompLabel & compID
							, const size_t label)  
	{
		assert( curNode < compID.size() /* initial node out of node index range */);
		if (compID[curNode] != UINT_MAX)	// nothing to do
			return;
		
		OutEdge_iterator curEdge = g.getOutEdgesBegin(curNode);
		  // temporary stack and initialization
		typedef std::pair<OutEdge_iterator,OutEdge_iterator> IT;
		std::stack< IT > path;
		compID[curNode] = label;
		  // add first node's edge iteration
		path.push( IT( curEdge, g.getOutEdgesEnd(curNode) ) );
		
		  // label component of current node
		while ( !path.empty() ) {
			IT t = path.top(); // get head element
			path.pop();
			  // process this neighborhood iteration
			for (; t.first != t.second; ++(t.first)) {
				  // add edge target to component if not already existing
				if ( compID[t.first->getToIndex()] == UINT_MAX) {
					compID[t.first->getToIndex()] = label;	// set label
					  // add adjacency of target node to processing queue
					path.push( IT(g.getOutEdgesBegin(t.first->getToIndex()), g.getOutEdgesEnd(t.first->getToIndex())) );
				}
			}
		}
	}


////////////////////////////////////////////////////////////////////////////

	Graph_Interface::Edge_iterator
	Graph_Interface::
	getEdgesBegin(const IndexType & i, const IndexType & j) const
	{
		  // create edge iterator based on out edge iterators
		return Edge_iterator( getOutEdgesBegin(i), getOutEdgesEnd(i), j );
	}

////////////////////////////////////////////////////////////////////////////

	Graph_Interface::Edge_iterator
	Graph_Interface::
	getEdgesEnd(const IndexType & i, const IndexType & j) const
	{
		  // create iterator end based on out edge iterators
		return Edge_iterator( getOutEdgesEnd(i), getOutEdgesEnd(i), j );
	}

////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	  // Default construction
	Graph_Interface::
	Edge_iterator::
	Edge_iterator()
	 :	curEdge()
		, edgeEnd()
		, toIndex(UINT_MAX)
	{
	}

	///////////////////////////////////////////////////////////////////////////

	  // Construction by data

	Graph_Interface::
	Edge_iterator::
	Edge_iterator( const OutEdge_iterator & curEdge_
				, const OutEdge_iterator & edgeEnd_
				, const IndexType & toIndex_ )
	 :	curEdge( curEdge_ )
		, edgeEnd( edgeEnd_ )
		, toIndex( toIndex_ )
	{
		  // go to first edge targeting the wanted index if existing
		while( curEdge != edgeEnd && curEdge->to != toIndex ) {
			++curEdge;
		}
	}


	///////////////////////////////////////////////////////////////////////////

	  // Increment to next element
	Graph_Interface::
	Edge_iterator &
	Graph_Interface::
	Edge_iterator::
	operator ++ ()
	{
		  // check if end reached
		if (curEdge == edgeEnd)
			return *this;

		  // else increment
		++curEdge;
		  // go to first edge targeting the wanted index if existing
		while( curEdge != edgeEnd && curEdge->to != toIndex ) {
			++curEdge;
		} // report changed status
		return *this;
	}

	///////////////////////////////////////////////////////////////////////////

	  // Increment to next element
	Graph_Interface::
	Edge_iterator
	Graph_Interface::
	Edge_iterator::
	operator ++ (int)
	{
		  // check if end reached
		if (curEdge == edgeEnd)
			return *this;

		  // else create temporary copy of current state to be returned at the end
		Edge_iterator oldStatus(*this);
		  // increment
		++curEdge;
		  // go to first edge targeting the wanted index if existing
		while( curEdge != edgeEnd && curEdge->to != toIndex ) {
			++curEdge;
		}
		  // access to the old non-altered iterator
		return oldStatus;
	}

	///////////////////////////////////////////////////////////////////////////

	  // access to the descriptive data of the current edge
	  // @return the current edge data
	const
	Graph_Interface::
	EdgeDescriptor&
	Graph_Interface::
	Edge_iterator::
	operator * () const
	{
		  // access to the dereferenced data object
		return *(curEdge);
	}

	///////////////////////////////////////////////////////////////////////////

	  // access to the descriptive data of the current edge
	  // @return the current edge data
	const
	Graph_Interface::
	EdgeDescriptor* const
	Graph_Interface::
	Edge_iterator::
	operator -> () const
	{
		  // access to the data object
		return curEdge.operator->();
	}

	///////////////////////////////////////////////////////////////////////////

	  // Equality comparison
	  // @return true if both iterators are pointing to the same edge
	bool
	Graph_Interface::
	Edge_iterator::
	operator == ( const Edge_iterator & it ) const
	{
		  // compare EdgeDescriptors
		return this->curEdge == it.curEdge;
	}

	///////////////////////////////////////////////////////////////////////////

	  // Inequality comparison
	  // @return true if both iterators are NOT pointing to the same edge
	bool
	Graph_Interface::
	Edge_iterator::
	operator != ( const Edge_iterator & it ) const
	{
		return this->curEdge != it.curEdge;
	}

	///////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////



} // namespace sgm

std::ostream&
operator <<( std::ostream & out, const sgm::Graph_Interface& g ) 
{
	for (size_t i=0; i<g.getNodeNumber(); ++i) {
		out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
		sgm::Graph_Interface::OutEdge_iterator eEnd = g.getOutEdgesEnd(i);
		for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i);
				e != eEnd; ++e)
		{
			out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
		}
		out <<" |\n";
	}
	return out;
}
