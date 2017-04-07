#ifndef GGL_RULEGRAPH_HH_
#define GGL_RULEGRAPH_HH_

	
////////////////////////////////////////////////////////////////////////////////

#include "ggl/Rule.hh"

#include <sgm/Graph_Interface.hh>
#include <sgm/Pattern.hh>
#include <sgm/PA_OrderCheck.hh>

#include <set>

namespace ggl {



	  /*! @brief Rule left side pattern
	   *
	   *  Implements the sgm::Pattern_Interface class to allow for the search of
	   *  the left side of a given ggl::Rule object using the sgm library.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class LeftSidePattern : public sgm::Graph_Interface, public sgm::Pattern_Interface {
		
	public:
		
		typedef std::set<IndexType> IndexSet;
		
	protected:

		  //! handle to the property map used to derive the edge label
		typedef boost::property_map<Rule::CoreGraph, Rule::EdgeLabelProperty>::const_type  EdgeLabelMap;
		  //! handle to the property map used to derive the node label
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeLabelProperty>::const_type  NodeLabelMap;
		  //! handle to the property map used to derive the local node indices
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeIndexProperty>::const_type  NodeIndexMap;
		  //! handle to the property map used to derive the local edge context
		typedef boost::property_map<Rule::CoreGraph, Rule::EdgeContextProperty>::const_type  EdgeContextMap;
		  //! handle to the property map used to derive the local node context
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeContextProperty>::const_type  NodeContextMap;

		  //! The Rule object to map into the Graph_Interface.
		const Rule& rule;
		
		  //! Holds the first index of each component of this graph
		  //! Has to be sorted in increasing indices !!!
		IndexSet firstOfEachComponent;
		
		  //! Compontent labeling
		sgm::Graph_Interface::CompLabel compLabels;
		
		  //! the list of order checks to apply to break all symmetries for this
		  //! rule 
		mutable sgm::PA_OrderCheck::CheckList* symmBreakList;
		
		  //! direct access of the edge label map of the core graph
		EdgeLabelMap edgeLabel;
		  //! direct access of the node label map of the core graph
		NodeLabelMap nodeLabel;
		  //! direct access of the node index map of the core graph
		NodeIndexMap nodeIndex;
		  //! direct access of the edge context map of the core graph
		EdgeContextMap edgeContext;
		  //! direct access of the node context map of the core graph
		NodeContextMap nodeContext;

	public:
		
		  //! Constructs a LeftSidePattern for the given rule
		  //! @param rule the Rule this object will be a left side pattern of
		LeftSidePattern( const Rule& rule );
		
		  //! Copy construction.
		  //! @param toCopy the  left side pattern to copy
		LeftSidePattern( const LeftSidePattern& toCopy );
		
		virtual
		~LeftSidePattern( );
		
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
		
		  //! Constant access to the internal Rule object.
		  //! @return the internal Rule object
		virtual
		const Rule &
		getRule(void) const;
		
		  /*! Creates a Pattern_Automorphism object to handle symmetries when
		   * matching the LeftSidePattern of this rule.
		   * The GRAPHMATCHER template argument is used to choose a
		   * sgm::GraphMatcher class to be used for automorphism detection.
		   * @return the Pattern_Automorphism for this rule
		   */
		template< class GRAPHMATCHER >
		sgm::PA_OrderCheck
		getGraphAutomorphismT( ) const;

		  /*! Creates a Pattern_Automorphism object to handle symmetries when
		   * matching the LeftSidePattern of this rule. Calls
		   * getGraphAutomorphismT< sgm::GM_vf2 >().
		   * @return the Pattern_Automorphism for this rule
		   */
		virtual
		sgm::PA_OrderCheck
		getGraphAutomorphism( ) const;
		
		  /*! Access to the additional constraints that have to be fulfilled by
		   * a valid matching of the left side of a graph grammar rule.
		   * @return the additional sgm::Pattern_Interface::Match_Constraint's to be matched
		   */
		virtual
		const ConstraintVec &
		getConstraints( void ) const;

		  /*! Access to *this, since this is the pattern graph..
		   * @return the pattern graph to be matched
		   */
		virtual
		const Graph_Interface &
		getPatternGraph( void ) const;

		  //! Access to the wildcard to be used when matching this pattern onto
		  //! some other graph.
		  //! @return the wildcard string to be used for edge and node labels,
		  //!         or NULL if no wildcard should be applied
		virtual
		const std::string *
		getUsedWildcard( void ) const;

		
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
			Rule::CoreGraph::out_edge_iterator cur_edge;
			  //! the edge iterator that marks the end of
			  //! the local adjacency list
			Rule::CoreGraph::out_edge_iterator list_end;

			  //! access to the parent instance for fast graph access
			const LeftSidePattern & parent;

		public:

			  //! Construction of a new edge descriptor given the our edge
			  //! iterators
			  //! @param cur_edge the current boost edge iterator
			  //! @param neigh_end the end of the edge iteration
			  //! @param parent access to the parent instance for graph access
			EdgeDescriptor(	const Rule::CoreGraph::out_edge_iterator& cur_edge,
							const Rule::CoreGraph::out_edge_iterator& neigh_end,
							const LeftSidePattern & parent );

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
			bool
			operator==(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
			bool
			operator!=(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
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

		  /*!
		   * Access to the component labeling of this pattern.
		   * @return the component labeling
		   */
		const CompLabel&
		getComponentLabeling( void ) const;
		
		  /*!
		   * Access to the first index of each component of this pattern.
		   * @return a set of the first index for each component
		   */
		const IndexSet&
		getFirstOfEachComponent( void ) const;
		
	};

	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ggl

	
////////////////////////////////////////////////////////////////////////////////

namespace ggl {


	  /*! @brief Rule right side graph
	   *
	   *  Implements the sgm::Graph_Interface class to allow for the search of
	   *  the right side of a given Rule object using the sgm library.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class RightSidePattern : public sgm::Graph_Interface {
		
	protected:
		
		  //! handle to the property map used to derive the edge label
		typedef boost::property_map<Rule::CoreGraph, Rule::EdgeLabelProperty>::const_type  EdgeLabelMap;
		  //! handle to the property map used to derive the node label
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeLabelProperty>::const_type  NodeLabelMap;
		  //! handle to the property map used to derive the right node label
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeRightLabelProperty>::const_type  NodeRightLabelMap;
		  //! handle to the property map used to derive the local node indices
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeIndexProperty>::const_type  NodeIndexMap;
		  //! handle to the property map used to derive the local edge context
		typedef boost::property_map<Rule::CoreGraph, Rule::EdgeContextProperty>::const_type  EdgeContextMap;
		  //! handle to the property map used to derive the local node context
		typedef boost::property_map<Rule::CoreGraph, Rule::NodeContextProperty>::const_type  NodeContextMap;


		  //! The Rule object to map into the Graph_Interface.
		const Rule& rule;
		

		  //! direct access of the edge label map of the core graph
		EdgeLabelMap edgeLabel;
		  //! direct access of the node label map of the core graph
		NodeLabelMap nodeLabel;
		  //! direct access of the node label map of the core graph
		NodeRightLabelMap nodeRightLabel;
		  //! direct access of the node index map of the core graph
		NodeIndexMap nodeIndex;
		  //! direct access of the edge context map of the core graph
		EdgeContextMap edgeContext;
		  //! direct access of the node context map of the core graph
		NodeContextMap nodeContext;

	public:
		
		RightSidePattern( const Rule& rule );
		
		virtual
		~RightSidePattern( );
		
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
		
		  //! Constant access to the internal Rule object.
		  //! @return the internal Rule object
		virtual
		const Rule &
		getRule(void) const;
		
		
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
			Rule::CoreGraph::out_edge_iterator cur_edge;
			  //! the edge iterator that marks the end of
			  //! the local adjacency list
			Rule::CoreGraph::out_edge_iterator list_end;

			  //! access to the parent instance for fast graph access
			const RightSidePattern & parent;

		public:

			  //! Construction of a new edge descriptor given the our edge
			  //! iterators
			  //! @param cur_edge the current boost edge iterator
			  //! @param neigh_end the end of the edge iteration
			  //! @param parent access to the parent instance for graph access
			EdgeDescriptor(	const Rule::CoreGraph::out_edge_iterator& cur_edge,
							const Rule::CoreGraph::out_edge_iterator& neigh_end,
							const RightSidePattern & parent );

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
			bool
			operator==(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
			bool
			operator!=(const EdgeDescriptor& ed ) const;

			  //! Inequality comparison
			  //! @param ed the edge to compare to
			  //! @return true if both descriptors describe different edges
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

		
	};

	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
	
} // namespace ggl

// template implementation
#include "ggl/RuleGraph.icc"


#endif /*RULEGRAPH_HH_*/
