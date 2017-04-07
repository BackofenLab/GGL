#ifndef SGM_GRAPH_BOOSTV_P_HH_
#define SGM_GRAPH_BOOSTV_P_HH_

#include "sgm/Graph_Interface.hh"
#include "sgm/Graph_Interface.hh"

#include <boost/graph/adjacency_list.hpp>

#include <vector>
#include <utility>

namespace sgm {


#define SGM_GRAPH_BOOSTV_P_TEMPLATE \
	template <	class GRAPH, \
				typename NODE_LABEL_PROPERTY, \
				typename EDGE_LABEL_PROPERTY, \
				typename NODE_INDEX_PROPERTY \
			>

#define SGM_GRAPH_BOOSTV_P_TYPE \
	Graph_boostV_p <	GRAPH \
						, NODE_LABEL_PROPERTY \
						, EDGE_LABEL_PROPERTY \
						, NODE_INDEX_PROPERTY \
					>

	  /*! @brief Graph_Interface wrapper for set of boost graphs
	   *
	   *  This class implements the Graph interface of a labeled undirected
	   *  graph. It is used as an interface of a vector of POINTERS of
	   *  undirected labeled
	   *  boost graphs to the sub graph matching algorithms of the library.
	   *  The template parameter GRAPH represents the type of the
	   *  boost graphs to use. The type names NODE_LABEL_PROPERTY and
	   *  EDGE_LABEL_PROPERTY are used to generate the node and edge labels of
	   *  the  graph. The type name NODE_INDEX_PROPERTY is used to derive the
	   *  corresponding local node index for a certain vertex in a graph.
	   *  The mapped set of graphs is represented as ONE unconnected
	   *  graph to the library.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	template <	class GRAPH, 
				typename NODE_LABEL_PROPERTY = boost::vertex_name_t, 
				typename EDGE_LABEL_PROPERTY = boost::edge_name_t,
				typename NODE_INDEX_PROPERTY = boost::vertex_index_t >
	class Graph_boostV_p : public Graph_Interface {
		
	protected:

		  //! index property map
		typedef typename boost::property_map<GRAPH, NODE_INDEX_PROPERTY>::const_type NodeIndexMap;

		  //! edge label property map
		typedef typename boost::property_map<GRAPH, EDGE_LABEL_PROPERTY>::const_type EdgeLabelMap;
		
		  //! the vector of graphs to map
		const std::vector< const GRAPH * > & graphs;
		
		  //! the sum of node numbers up to graphs[i]
		std::vector< size_t > graphSize;
		
	public:
		
		  //! the location information of a node in the global index to the
		  //! graph index in graphs and the local node index in that graph.
		  //! LocalIndex.first gives the graph index in [0,graphs.size())
		  //! LocalIndex.second gives the local node in [0,graphs[LI.first].nodes)
		typedef std::pair< size_t, size_t > LocalIndex;
		
		typedef GRAPH InternalBoostGraph;
		
		
		  //! Construction of the interface graph.
		  //! @param graphs the vector of boost graph pointers to use as 
		  //!               internal data structure (all != NULL)
		Graph_boostV_p( const std::vector< const GRAPH * > & graphs );
		
		  //! Destruction
		virtual 
		~Graph_boostV_p();
		
		  //! Access to the number of nodes of the graph
		  //! @return the overall node number 
		virtual
		size_t
		getNodeNumber(void) const;
		
		  //! Access to the local index in graphs of the virtual global node i
		  //! @param i index of the global node of intrest
		  //! @return the local index of i in graphs and the local index of i in
		  //!         that graph
		LocalIndex
		getLocalIndex( const IndexType & i ) const;
		
		  //! Access to the label of a specified node
		  //! @param i the index of the node of interest
		  //! @return a string representation of the node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const;
		
		
		  //! Access to the internal list of GRAPH object pointers
		  //! @return a reference to the graphs object
		const std::vector< const GRAPH* > & 
		getGraphs(void) const;


	public:
		
		  /*!
		   * Special edge descriptor to enable the iteration over the edges
		   * contained in the internally used boost graph instance.
		   *
		   */
		class EdgeDescriptor : public Graph_Interface::EdgeDescriptor {
		protected:

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

			  //! access to the graph instance
			const GRAPH& graph;

			  //! the index shift to correct node indices to the global index
			const size_t nodeIndexShift;

			  //! node index access
			NodeIndexMap nodeIndexMap;
			  //! edge label access
			EdgeLabelMap edgeLabelMap;

		public:

			  //! Construction of a new edge descriptor given the our edge
			  //! iterators
			  //! @param cur_edge the current boost edge iterator
			  //! @param neigh_end the end of the edge iteration
			  //! @param graph the graph instance the edge iterator is based on
			  //! @param nodeIndexShift index shift needed to correct node
			  //!        indices to the global index
			EdgeDescriptor(	const typename GRAPH::out_edge_iterator& cur_edge,
							const typename GRAPH::out_edge_iterator& neigh_end,
							const GRAPH& graph,
							const size_t nodeIndexShift
							);

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
		  //! in the graphs. For each graph it opens the property map for
		  //! NODE_INDEX_PROPERTY and sets an index with [0,num_vertices(g))
		  //! along the ordering defined by the vertex iterators.
		  //! Therefore, this function can be used to set the necessary local 
		  //! indexes before creating a Graph_boostV_p object.
		  //! NOTE : The property map is NOT CONSTRUCTED has to be already part
		  //! of the GRAPH definition! It is ONLY FILLED!
		static 
		void
		indexGraphs ( std::vector< GRAPH * > & graphs );
	};

} // namespace sgm


  // include template implementation of the class
#include "sgm/Graph_boostV_p.icc"


#endif /*GRAPH_BOOSTV_P_HH_*/
