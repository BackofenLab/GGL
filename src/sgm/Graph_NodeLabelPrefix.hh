#ifndef SGM_GRAPH_SUBNODELABEL_HH_
#define SGM_GRAPH_SUBNODELABEL_HH_

#include "sgm/Graph_Interface.hh"

namespace sgm {

	  /*! @brief Wrapper for a graph with reduced node label information
	   *
	   *  For all node labels, only the prefix up to a given separator is
	   *  returned as new node label
	   *
	   *  @author Martin Mann (c) 2017 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Graph_NodeLabelPrefix : public sgm::Graph_Interface {

	protected:
		
		//! wrapped fullGraph representation to reduce the node label information
		const sgm::Graph_Interface & fullGraph;
		
		//! the regular expression pattern to be applied to all node labels
		const std::string nodeLabelSeparator;

	public:

		//! Construction
		//! @param fullGraph the graph representing to wrap;
		//!    node labels are reduced to their concatenation of matched
		//!    group substrings for the given node label pattern
		//! @param nodeLabelSeparator the (non-empty) separator that marks the beginning
		//!    of the node label part to be ignored; everything IN FRONT of
		//!    this substring is returned as node label;
		Graph_NodeLabelPrefix( const sgm::Graph_Interface & fullGraph,
				const std::string & nodeLabelSeparator )
		 :	fullGraph(fullGraph)
			, nodeLabelSeparator(nodeLabelSeparator)
		{
			if (nodeLabelSeparator.size()==0) {
				throw std::runtime_error("Graph_NodeLabelPrefix() : nodeLabelSeparator empty");
			}
		}

		  //! Destruction
		virtual 
		~Graph_NodeLabelPrefix() {}
		
		  //! Access to the number of nodes of the graph
		  //! @return the overall node number 
		virtual
		size_t
		getNodeNumber(void) const
		{ return fullGraph.getNodeNumber(); }
		
		  //! Access to iteration begin for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator to the first edge within the adjacency of i
		virtual
		OutEdge_iterator
		getOutEdgesBegin( const IndexType & i ) const
		{ return fullGraph.getOutEdgesBegin( i ); }

		  //! Access to iteration end for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator the end of the adjacency iteration of i
		virtual
		OutEdge_iterator
		getOutEdgesEnd( const IndexType & i ) const
		{ return fullGraph.getOutEdgesEnd( i ); }
		
		  //! Access to the label of a specified atom node where the original
		  //!    atom label is reduced to their substrings up to the first
		  //!    occurrence of the class label separator ':'
		  //! @param i the index of the node of interest
		  //! @return a string representation of the atom node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const
		{
			std::string origLabel = fullGraph.getNodeLabel(i);
			// return prefix of original atom label UP TO class label separator ':'
			return origLabel.substr( 0, origLabel.find(nodeLabelSeparator) );
		}

		//! Provides access to the original graph without reduced atom labels.
		//! @return the original graph wrapped by this graph
		const sgm::Graph_Interface &
		getWithFullAtomLabels() const
		{ return fullGraph; }

	};

} // namespace sgm



#endif /*SGM_GRAPH_SUBNODELABEL_HH_*/
