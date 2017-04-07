#ifndef SGM_GRAPH_INTERFACE_HH_
#define SGM_GRAPH_INTERFACE_HH_

#include <vector>
#include <string>
#include <memory>
#include <climits>

namespace sgm {

	  /*! @brief Interface of graphs for graph matching
	   *
	   *  The class is used as an abstract interface to allow for a generic
	   *  interface of the search algorithms in the library.
	   *  It describes an undirected labeled graph.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Graph_Interface {
	public:

		  //! Type for global indicex of nodes in the graph.
		typedef size_t IndexType;
		
		class Edge_iterator;
		class OutEdge_iterator;

		  //! Container that stores the global edge information, i.e. from, to,
		  //! and the label of the edge.
		class EdgeDescriptor {
		public:
			  //! Default Construction
			EdgeDescriptor();

			  //! Construction
			  //! @param from the source of the edge
			  //! @param to the target of the edge
			EdgeDescriptor(	const Graph_Interface::IndexType& from, 
							const Graph_Interface::IndexType& to );

			  //! Destruction
			virtual
			~EdgeDescriptor();

		protected:

			  // define friends for direct data member access
			friend class Edge_iterator;
			  // define friends for direct data member access
			friend class OutEdge_iterator;

			  //! the source of the edge described
			Graph_Interface::IndexType from;
			  //! the target of the edge described
			Graph_Interface::IndexType to;
			
		public:
			  //! Access to the source of the edge.
			  //! @return the global source node index
			const Graph_Interface::IndexType&
			getFromIndex(void) const;
			  //! Access to the target of the edge.
			  //! @return the global target node index
			const Graph_Interface::IndexType&
			getToIndex(void) const;
			  //! Access to the label of the edge.
			  //! @return the edge label
			virtual
			const std::string&
			getEdgeLabel(void) const = 0;
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
			operator!=(const EdgeDescriptor& ed ) const = 0;

			  //! Iterator support
			  //! @return the next EdgeDescriptor in the adjacency
			virtual
			EdgeDescriptor&
			operator++() = 0;

			  //! Create a heap copy of this object. NOTE: this has to be
			  //! removed by the calling function.
			  //! @return a new heap copy of this object
			virtual
			EdgeDescriptor *
			clone() const = 0;
		};
		
		  /*!
		   * Generic constant iterator class to enumerate the outgoing edges
		   * for a given node. For each out edge an EdgeDescriptor instance is
		   * provided.
		   */
		class OutEdge_iterator {
		protected:

			  //! the data object to be used for iteration. Note: the object
			  //! can be also any subclass of EdgeDescriptor.
			std::auto_ptr< EdgeDescriptor > data;

		public:

			  //! Default construction
			OutEdge_iterator();

			  //! Construction from data
			  //! @param data the EdgeDescriptor data to use
			OutEdge_iterator( const EdgeDescriptor & data );

			  //! Copy construction
		      //! @param toCopy the object to make this iterator a copy of
			OutEdge_iterator( const OutEdge_iterator & toCopy );

			  //! Destruction
			~OutEdge_iterator();

			  //! Assign to the given iterator
		      //! @param toCopy the object to make this iterator a copy of
			OutEdge_iterator & operator = (const OutEdge_iterator & toCopy);

			  //! Prefix increment (++I) to next element
			OutEdge_iterator & operator ++ ();

			  //! Postfix increment (I++) to next element
			OutEdge_iterator operator ++ (int);

			  //! access to the descriptive data of the current edge
		      //! @return the current edge data
			const EdgeDescriptor& operator * () const;

			  //! access to the descriptive data of the current edge
		      //! @return the current edge data
			const EdgeDescriptor* const operator -> () const;

			  //! Equality comparison
		      //! @return true if both iterators are pointing to the same edge
			bool operator == ( const OutEdge_iterator & it ) const ;

			  //! Inequality comparison
		      //! @return true if both iterators are NOT pointing to the same edge
			bool operator != ( const OutEdge_iterator & it ) const ;

		};


		  /*!
		   * Generic constant iterator class to enumerate the edges between two
		   * nodes. It uses the OutEdge_iterator interface to perform the
		   * enumeration. For each edge an EdgeDescriptor is provided.
		   */
		class Edge_iterator {
		protected:

			  //! The OutEdge_iterator used to enumerate the edges
			OutEdge_iterator curEdge;
			  //! The OutEdge_iterator marking the end of the enumeration
			OutEdge_iterator edgeEnd;
			  //! The targeted index of the edges to be reported
			const size_t toIndex;

		public:

			  //! Default construction
			Edge_iterator();

			  //! Construction from data
			  //! @param curEdge the begin of the outEdges to iterator
			  //! @param edgeEnd the end of the outEdges to iterator
			  //! @param toIndex the targeted index of the edges to be reported
			Edge_iterator(	const OutEdge_iterator & curEdge
							, const OutEdge_iterator & edgeEnd
							, const IndexType & toIndex );

			  //! Prefix increment (++I) to next element
			Edge_iterator & operator ++ ();

			  //! Postfix increment (I++) to next element
			Edge_iterator operator ++ (int);

			  //! access to the descriptive data of the current edge
		      //! @return the current edge data
			const EdgeDescriptor& operator * () const;

			  //! access to the descriptive data of the current edge
		      //! @return the current edge data
			const EdgeDescriptor* const operator -> () const;

			  //! Equality comparison
		      //! @return true if both iterators are pointing to the same edge
			bool operator == ( const Edge_iterator & it ) const ;

			  //! Inequality comparison
		      //! @return true if both iterators are NOT pointing to the same edge
			bool operator != ( const Edge_iterator & it ) const ;

		};

		  //! Destruction
		virtual 
		~Graph_Interface() {}
		
		  //! Access to the number of nodes of the graph
		  //! @return the overall node number 
		virtual
		size_t
		getNodeNumber(void) const = 0;
		
		  //! Access to iteration begin for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator to the first edge within the adjacency of i
		virtual
		OutEdge_iterator
		getOutEdgesBegin( const IndexType & i ) const = 0;

		  //! Access to iteration end for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator the end of the adjacency iteration of i
		virtual
		OutEdge_iterator
		getOutEdgesEnd( const IndexType & i ) const = 0;
		
		  //! Access to the label of a specified node
		  //! @param i the index of the node of interest
		  //! @return a string representation of the node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const = 0;
		
		  //! Access to the label of a specified edge
		  //! @param i the index of the first end node of interest
		  //! @param j the index of the second end node of interest
		  //! @return the edge iterator pointing to the first edge between
		  //! node i and j or to getEdgesEnd(i,j) if no edge exists
		virtual
		Edge_iterator
		getEdgesBegin(const IndexType & i, const IndexType & j) const;

		  //! Access to the label of a specified edge
		  //! @param i the index of the first end node of interest
		  //! @param j the index of the second end node of interest
		  //! @return the edge iterator pointing to the first edge between
		  //! node i and j or to getEdgesEnd(i,j) if no edge exists
		virtual
		Edge_iterator
		getEdgesEnd(const IndexType & i, const IndexType & j) const;

		  //! Interface equality comparison : NOTE : this function checks
		  //! whether the interfaces of two graphs are identical or not, i.e.
		  //! nodes with same index are equal and the order of adjacent edges
		  //! is the same. This function performs NO GRAPH ISOMORPHISM !!!
		  //! Thus, slight changes in the node order etc. will result in a
		  //! non-equal interface!
		  //! @param toCompare the Graph_Interface to compare to
		  //! @return true if both graph interfaces are equal
		virtual
		bool 
		operator==(const Graph_Interface& toCompare ) const;

		  //! Interface inequality comparison : NOTE : this function checks
		  //! whether the interfaces of two graphs are identical or not, i.e.
		  //! nodes with same index are equal and the order of adjacent edges
		  //! is the same. This function performs NO GRAPH ISOMORPHISM !!!
		  //! Thus, slight changes in the node order etc. will result in a
		  //! non-equal interface!
		  //! @param toCompare the Graph_Interface to compare to
		  //! @return true if both graph interfaces are different
		virtual
		bool 
		operator!=(const Graph_Interface& toCompare ) const;

		//////////////////  STATIC MEMBERS /////////////////////////////////////
		
		typedef std::vector<size_t> CompLabel;
		
		  /*!
		   * Computes the connected components of a given graph. It stores the
		   * connected component ID of each node within the provided container
		   * and returns the number of components.
		   * @param g the graph to check for connected components
		   * @param compID the container to store the component ID information
		   * within
		   * @return the number of connected components within g
		   */
		static
		size_t
		connectedComponents(	const Graph_Interface& g
								, CompLabel & compID );
		
	protected:
		

		  /*!
		   * Performs a depths-first-search labeling of all accessible nodes
		   * starting from curNode within g. The reached nodes are labeled
		   * with the given label in the component ID container compID
		   *
		   * @param the graph to screen
		   * @param curNode the node to start the search from
		   * @param compID the container to store the component ID information
		   * within
		   * @param label the component ID to be used for all nodes of the
		   * connected component reachable from curNode
		   */
		static
		void 
		labelAdjacentNodes( const Graph_Interface& g
							, const size_t curNode
							, CompLabel & compID
							, const size_t label );  

	};

} // namespace sgm


#include <iostream>
#include <iomanip>


  /*! Prints a Graph_Interface instance to stream. For each node its label and
   * the adjancent nodes including the edge label is printed.
   * 
   * @param out the stream to write to
   * @param g the graph to write
   * @return the modified out stream
   */
std::ostream&
operator <<( std::ostream & out, const sgm::Graph_Interface& g );


 // include inline function bodies
#include "sgm/Graph_Interface.icc"

#endif /*SGM_GRAPH_INTERFACE_HH_*/
