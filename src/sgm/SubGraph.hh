#ifndef SGM_SUBGRAPH_HH_
#define SGM_SUBGRAPH_HH_

#include "sgm/Graph_Interface.hh"
#include "sgm/Match.hh"

namespace sgm
{
	
	/*! @brief Subgraph of a Graph_Interface for matching
	 *
	 * Represents only a subset of the nodes of a graph and all edges between
	 * these nodes by wrapping the source graph.
	 * 
	 * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	class SubGraph : public sgm::Graph_Interface
	{
	public:
		
		  //! Type that defines a node index subset to be present in a subgraph
		typedef std::vector< IndexType > NodeList;
		
	protected:

		  /*!
		   * Special edge descriptor to enable the iteration over the edges
		   * covered by this subgraph.
		   *
		   */
		class EdgeDescriptor : public sgm::Graph_Interface::EdgeDescriptor {
		protected:

			  //! the source of the edge described
			using Graph_Interface::EdgeDescriptor::from;
			  //! the target of the edge described
			using Graph_Interface::EdgeDescriptor::to;

			  //! the iterator to the current edge
			OutEdge_iterator curEdge;
			  //! the iterator to the edge iteration end
			OutEdge_iterator edgeEnd;
			  //! access to the parent class for full subgraph information
			const SubGraph & parent;

		public:
			  //! Construction
			  //! @param curEdge the iterator to the current edge
			  //! @param edgeEnd the iterator to the iteration end
			  //! @param parent access to the parent full graph information
			EdgeDescriptor(	const Graph_Interface::OutEdge_iterator & curEdge
							, const Graph_Interface::OutEdge_iterator & edgeEnd
							, const SubGraph & parent );

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
		

	protected:
		
		static const IndexType NOT_MAPPED_INDEX;
		
		  //! The original graph this object is a subgraph of.
		const Graph_Interface & fullGraph;
		
		  //! The subset of node indices of fullGraph represented by this graph
		NodeList nodeList;
		
		  //! Mapping of global indices of fullGraph to local indices of this
		  //! subgraph
		Match full2sub;
		
	public:
		
		  /*!
		   * Construction of a subgraph from explicit node list.
		   * 
		   * @param fullGraph the original graph this object should be a 
		   *                  subgraph of
		   * @param nodeList the subset of node indices present in the subgraph
		   */
		SubGraph(	const Graph_Interface & fullGraph
					, const NodeList & nodeList );
		
		  /*!
		   * Construction of a subgraph from a component labeling.
		   * 
		   * @param fullGraph the original graph this object should be a 
		   *                  subgraph of
		   * @param compLabel component labeling of the fullGraph
		   * @param subGraphLabel the label in compLabel that defines the nodes
		   *                  of this subgraph
		   */
		SubGraph(	const Graph_Interface & fullGraph
					, const CompLabel & compLabel
					, const size_t subGraphLabel );
		
		
		  //! Destruction of the subgraph
		virtual ~SubGraph();

	
		
		  //! Access to the number of nodes of the subgraph
		  //! @return the overall node number 
		virtual
		size_t
		getNodeNumber(void) const;
		
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
		
		  //! Access to the label of a specified node
		  //! @param i the index of the node of interest
		  //! @return a string representation of the node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const;
		
};

} // namespace sgm

#include "sgm/Pattern.hh"

namespace sgm
{

	/*!
	 * A Pattern that represents only a subset of the nodes of a graph and all
	 * edges between these nodes by wrapping the source graph and some
	 * additional constraints
	 *
	 * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	class SubGraphPattern : public SubGraph, public Pattern
	{
	protected:

		  //! the graph to be represented as a pattern
		using Pattern::graph;

		  //! the additional match constraints to be fulfilled by each match
		using Pattern::matchConstraints;

		  //! the wildcard string to be used for matching
		using Pattern::usedWildcard;

	public:

		  /*!
		   * Construction of a subgraph from explicit node list.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param nodeList the subset of node indices present in the subgraph
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const NodeList & nodeList );

		  /*!
		   * Construction of a subgraph from a component labeling.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param compLabel component labeling of the fullGraph
		   * @param subGraphLabel the label in compLabel that defines the nodes
		   *                  of this subgraph
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const CompLabel & compLabel
							, const size_t subGraphLabel );


		  /*!
		   * Construction of a subgraph from explicit node list.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param nodeList the subset of node indices present in the subgraph
		   * @param matchConstraints the additional matching constraints to be
		   *        fulfilled for this pattern, NOTE: indices have to correspond
		   *        to the initial fullGraph instance and to be compatible with
		   *        the subgraph definition !
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const NodeList & nodeList
							, const ConstraintVec & matchConstraints );

		  /*!
		   * Construction of a subgraph from a component labeling.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param compLabel component labeling of the fullGraph
		   * @param subGraphLabel the label in compLabel that defines the nodes
		   *                  of this subgraph
		   * @param matchConstraints the additional matching constraints to be
		   *        fulfilled for this pattern, NOTE: indices have to correspond
		   *        to the initial fullGraph instance and to be compatible with
		   *        the subgraph definition !
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const CompLabel & compLabel
							, const size_t subGraphLabel
							, const ConstraintVec & matchConstraints );

		  /*!
		   * Construction of a subgraph from explicit node list.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param nodeList the subset of node indices present in the subgraph
		   *  @param wildcardToUse the wildcard to use for matching
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const NodeList & nodeList
							, const std::string &  wildcardToUse );

		  /*!
		   * Construction of a subgraph from a component labeling.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param compLabel component labeling of the fullGraph
		   * @param subGraphLabel the label in compLabel that defines the nodes
		   *                  of this subgraph
		   *  @param wildcardToUse the wildcard to use for matching
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const CompLabel & compLabel
							, const size_t subGraphLabel
							, const std::string &  wildcardToUse );


		  /*!
		   * Construction of a subgraph from explicit node list.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param nodeList the subset of node indices present in the subgraph
		   * @param matchConstraints the additional matching constraints to be
		   *        fulfilled for this pattern, NOTE: indices have to correspond
		   *        to the initial fullGraph instance and to be compatible with
		   *        the subgraph definition !
		   *  @param wildcardToUse the wildcard to use for matching
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const NodeList & nodeList
							, const ConstraintVec & matchConstraints
							, const std::string &  wildcardToUse );

		  /*!
		   * Construction of a subgraph from a component labeling.
		   *
		   * @param fullGraph the original graph this object should be a
		   *                  subgraph of
		   * @param compLabel component labeling of the fullGraph
		   * @param subGraphLabel the label in compLabel that defines the nodes
		   *                  of this subgraph
		   * @param matchConstraints the additional matching constraints to be
		   *        fulfilled for this pattern, NOTE: indices have to correspond
		   *        to the initial fullGraph instance and to be compatible with
		   *        the subgraph definition !
		   *  @param wildcardToUse the wildcard to use for matching
		   */
		SubGraphPattern(	const Graph_Interface & fullGraph
							, const CompLabel & compLabel
							, const size_t subGraphLabel
							, const ConstraintVec & matchConstraints
							, const std::string &  wildcardToUse );


		  //! Destruction of the subgraph
		virtual ~SubGraphPattern();


	protected:

		 /*!
		  * Does the reindexing of the nodes within the given match constraints
		  * onto the new indices within the subgraph, where the index mapping is
		  * given by full2sub.
		  *
		  * @param originalConstraints the matching constraints to be remapped
		  * @param full2sub the mapping of the indices used in the given
		  *        constraints onto nodes within the subgraph.
		  * @param remappedConstraints the container that will hold the
		  *        remapped matching constraints
		  */
		static
		void
		remapConstraints( const ConstraintVec & originalConstraints
							, const Match & full2sub
							, ConstraintVec & remappedConstraints );

	};

} // namespace sgm

// inline implementations
#include "sgm/SubGraph.icc"

#endif /*SUBGRAPH_HH_*/
