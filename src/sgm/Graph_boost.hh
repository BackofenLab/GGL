#ifndef SGM_GRAPH_BOOST_HH_
#define SGM_GRAPH_BOOST_HH_

#include "sgm/Graph_Interface.hh"

#include <boost/graph/adjacency_list.hpp>

#include <vector>
#include <utility>


namespace sgm {


#define SGM_GRAPH_BOOST_TEMPLATE \
	template <	class GRAPH, \
				typename NODE_LABEL_PROPERTY, \
				typename EDGE_LABEL_PROPERTY, \
				typename NODE_INDEX_PROPERTY \
			>

#define SGM_GRAPH_BOOST_TYPE \
	sgm::Graph_boost <	GRAPH \
						, NODE_LABEL_PROPERTY \
						, EDGE_LABEL_PROPERTY \
						, NODE_INDEX_PROPERTY \
					>


	  /*! @brief Graph_Interface wrapper for boost graphs
	   *
	   *  This class implements the Graph interface of a labeled undirected
	   *  graph. It is used as an interface of an undirected labeled
	   *  boost graph to the sub graph matching algorithms of the sgm library.
	   *  The template parameter GRAPH represents the type of the
	   *  boost graph to use. The type names NODE_LABEL_PROPERTY and
	   *  EDGE_LABEL_PROPERTY are used to generate the node and edge labels of
	   *  the  graph. The type name NODE_INDEX_PROPERTY is used to derive the
	   *  corresponding local node index for a certain vertex in covered graph.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	template <	class GRAPH, 
				typename NODE_LABEL_PROPERTY = boost::vertex_name_t, 
				typename EDGE_LABEL_PROPERTY = boost::edge_name_t,
				typename NODE_INDEX_PROPERTY = boost::vertex_index_t >
	class Graph_boost : public Graph_Interface {
		
	protected:
		
		  //! the graph to map
		const GRAPH & graph;
		  //! the number of nodes of graph
		const size_t graphSize;
		
		typedef typename boost::property_map<GRAPH, NODE_INDEX_PROPERTY>::const_type NodeIndexMap;
		  //! node index access
		NodeIndexMap nodeIndexMap;

		  //! node label access
		typename boost::property_map<GRAPH, NODE_LABEL_PROPERTY>::const_type
			nodeLabelMap;

		typedef typename boost::property_map<GRAPH, EDGE_LABEL_PROPERTY>::const_type EdgeLabelMap;
		  //! edge label access
		EdgeLabelMap edgeLabelMap;

	public:
		
		typedef GRAPH InternalBoostGraph;
		
		  //! Construction of the interface graph.
		  //! @param graph the boost graph to use as internal data structure
		Graph_boost( const GRAPH & graph );
		
		  //! Destruction
		virtual 
		~Graph_boost();
		
		  //! Access to the number of nodes of the graph
		  //! @return the overall node number 
		virtual
		size_t
		getNodeNumber(void) const;
		
		  //! Access to the label of a specified node
		  //! @param i the index of the node of interest
		  //! @return a string representation of the node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const;
		
		  //! Access to the internal GRAPH object
		  //! @return a refence to the graph object
		const GRAPH & 
		getGraph(void) const;

	public:
		
		  /*!
		   * Special edge descriptor to enable the iteration over the edges
		   * contained in the internally used boost graph instance.
		   *
		   */
		class EdgeDescriptor : public Graph_Interface::EdgeDescriptor {
		protected:
			  //! the source of the edge described
			using Graph_Interface::EdgeDescriptor::from;
			  //! the target of the edge described
			using Graph_Interface::EdgeDescriptor::to;
			  //! the out edge iterator of the current target node within
			  //! the local adjacency list
			typename GRAPH::out_edge_iterator cur_edge;
			  //! the edge iterator that marks the end of
			  //! the local adjacency list
			typename GRAPH::out_edge_iterator list_end;

			  //! access to the parent instance for fast graph access
			const Graph_boost& parent;

		public:

			  //! Construction of a new edge descriptor given the our edge
			  //! iterators
			  //! @param cur_edge the current boost edge iterator
			  //! @param neigh_end the end of the edge iteration
			  //! @param parent access to the parent instance for graph access
			EdgeDescriptor(	const typename GRAPH::out_edge_iterator& cur_edge,
							const typename GRAPH::out_edge_iterator& neigh_end,
							const Graph_boost& parent );

			  //! Destruction
			virtual
			~EdgeDescriptor();

			  //! Access to the label of the edge.
			  //! @return the edge label
			virtual
			const std::string&
			getEdgeLabel(void) const;

			  //! Equality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe the same edge
			virtual
			bool
			operator==(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
			virtual
			bool
			operator!=(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
			virtual
			bool
			operator!=(const Graph_Interface::EdgeDescriptor& ed ) const;

			  //! Iterator support
			  //! @return the next EdgeDescriptor in the adjacency
			virtual
			EdgeDescriptor&
			operator++();

			  //! Create a heap copy of this object. NOTE: this has to be
			  //! removed by the calling function.
			  //! @return a new heap copy of this object
			virtual
			EdgeDescriptor *
			clone() const;

		};

		  //! Access to iteration begin for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator to the first edge within the adjacency of i
		virtual
		OutEdge_iterator
		getOutEdgesBegin( const IndexType & i ) const;

		  //! Access to iteration end for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator the end of the adjacency iteration of i
		virtual
		OutEdge_iterator
		getOutEdgesEnd( const IndexType & i ) const;

		//////////////////// STATIC MEMBERS ///////////////////////////////////

	public:
		
		  //! Static member function that allows for the indexing of the nodes
		  //! in the graph. It opens the property map for
		  //! NODE_INDEX_PROPERTY and sets an index with [0,num_vertices(g))
		  //! along the ordering defined by the vertex iterators.
		  //! Therefore, this function can be used to set the necessary local 
		  //! indexes before creating a Graph_boost object.
		  //! NOTE : The property map is NOT CONSTRUCTED has to be already part
		  //! of the GRAPH definition! It is ONLY FILLED!
		static 
		void
		indexGraph ( GRAPH & graph );

		  //! Counts the number of nodes a boost graph represents since the
		  //! number returned by boost::num_vertices might not be correct, e.g.
		  //! for boost::filtered_graph instances.
		  //! @param graph the graph to check
		  //! @return the graph's number of nodes
		static
		size_t
		countRealNodeNumber( const GRAPH & graph );

	};


} // namespace sgm


  // include template implementation of the class
#include "sgm/Graph_boost.icc"


#endif /*SGM_GRAPH_BOOST_HH_*/
