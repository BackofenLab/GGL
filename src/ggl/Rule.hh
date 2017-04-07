#ifndef GGL_RULE_HH_
#define GGL_RULE_HH_

#include "sgm/HashMap.hh"
#if HAVE_UNORDERED_MAP > 0
	#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP > 0
	#include <tr1/unordered_map>
#elif HAVE_GNU_HASH_MAP > 0
	#include <ext/hash_map>
#else
	#include <map>
#endif

#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "sgm/MC_Node.hh"
#include "sgm/MC_Edge.hh"

namespace ggl {

	

////////////////////////////////////////////////////////////////////////////////

	  /*! @brief Rule consistency codes
	   *
	   * A helper class that holds constants for consistency or error encodings
	   * that are independent of any template parameter.
	   *
	   * @author Martin Mann (c) 2012 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Rule_ConCodes {
	public:
		  //! consistency code : everything fine
		static const size_t C_Consistent;
		  //! consistency code : no rule ID is present
		static const size_t C_NoID;
		  //! consistency code : at least one edge shows wrong context usage
		static const size_t C_WrongEdgeContext;

	};

	
////////////////////////////////////////////////////////////////////////////////
	
	
	  /*! Graph grammar rule data
	   *
	   * A description of a graph grammar rule with left and righ rule side and
	   * context description.
	   *
	   * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */
	class Rule : public Rule_ConCodes {

	public:
		

		  //! Vector of node indices.
		typedef std::vector< size_t > NodeVec;

		  //! If a node has outgoing edges, this class contains the indices of the
		  //! ends of these edges.
#if HAVE_UNORDERED_MAP > 0
		typedef std::unordered_map< size_t, NodeVec > OutEdgeMap;
		typedef std::unordered_map< size_t, size_t > NodeMap;
#elif HAVE_TR1_UNORDERED_MAP > 0
		typedef std::tr1::unordered_map< size_t, NodeVec > OutEdgeMap;
		typedef std::tr1::unordered_map< size_t, size_t > NodeMap;
#elif HAVE_GNU_HASH_MAP > 0
		typedef __gnu_cxx::hash_map< size_t, NodeVec > OutEdgeMap;
		typedef __gnu_cxx::hash_map< size_t, size_t > NodeMap;
#else
		typedef std::map< size_t, NodeVec > OutEdgeMap;
		typedef std::map< size_t, size_t > NodeMap;
#endif

		  //! A type definition to encode the different use cases of a node or edge
		  //! of a graph grammar rule.
		enum RuleContext {  RULE_LEFT_SIDE		= 2 //!< graph element is ONLY part of the left rule side
							, RULE_RIGHT_SIDE	= 3 //!< graph element is ONLY part of the right rule side
							, RULE_LABEL_CHANGE	= (RULE_LEFT_SIDE * RULE_RIGHT_SIDE) //!< node is part of the both rule sides but the LABEL is CHANGING. NOTE: only valid for nodes!!!
							, RULE_CONTEXT		= 5 * RULE_LABEL_CHANGE  //!< graph element is part of the both rule sides = context of the rule
						};
		
		  //! This boost graph property is used to determine the index of a given
		  //! node along the iterator order.
		typedef boost::vertex_index_t	NodeIndexProperty;
		
		  //! This boost graph property is used to determine the label of a given
		  //! node. If the specified node changes the label in the left and right
		  //! side of the graph grammar rule, this property specifies the label of
		  //! the LEFT side of the node.
		typedef boost::vertex_name_t	NodeLabelProperty;

		  //! This boost graph property is used to determine if a node is part of
		  //! the left, right, or both sides of a graph grammar rule.
		struct NodeContextProperty {
		    typedef boost::vertex_property_tag kind;
		  };

		  //! This boost graph property is used to determine if a node is part of
		  //! the left, right, or both sides of a graph grammar rule.
		struct NodeRightLabelProperty {
		    typedef boost::vertex_property_tag kind;
		  };

		
		  //! The properties available for the nodes of a CoreGraph
		typedef	boost::property<	NodeIndexProperty, size_t
				, boost::property<	NodeLabelProperty, std::string
				, boost::property<	NodeRightLabelProperty, std::string
				, boost::property<	NodeContextProperty, RuleContext
					> > > >
						CoreGraph_NodeProperties;

		  //! This boost graph property is used to determine the label of a given
		  //! edge. If the specified edge changes the label in the left and right
		  //! side of the graph grammar rule, this property specifies the label of
		  //! the LEFT side of the edge.
		typedef boost::edge_name_t		EdgeLabelProperty;

		  //! This boost graph property is used to determine if an edge is part of
		  //! the left, right, or both sides of a graph grammar rule.
		struct EdgeContextProperty {
		    typedef boost::edge_property_tag kind;
		  };

		  /*! The properties available for the edges of a CoreGraph
		   *
		   * NOTE: EdgeContextProperty is NOT TO BE allowed to be RULE_LABEL_CHANGE.
		   * Label changes of edges are to be encoded by insertion and deletion of
		   * the edge, ie. present two times with according RULE_LEFT_SIDE and
		   * RULE_RIGHT_SIDE EdgeContextProperty value.
		   */
		typedef	boost::property<	EdgeLabelProperty, std::string
				, boost::property<	EdgeContextProperty, RuleContext
					> >
						CoreGraph_EdgeProperties;

		  //! The definition of a the internal graph representation of a graph
		  //! grammar rule of undirected graphs.
		typedef boost::adjacency_list<
		  boost::vecS,      				// store edges
		  boost::vecS,       				// store vertices
		  boost::undirectedS,				// is an undirected graph
		  CoreGraph_NodeProperties,     // (atom symbols etc)
		  CoreGraph_EdgeProperties      // (edge symbols etc)
		  >
			CoreGraph;

		  //! vector of constraints to be matched
		typedef std::vector< sgm::Pattern_Interface::Match_Constraint* >
			MatchConstraintVec;

	protected:
		/////////////////  DATA MEMBERS  /////////////////////////

		  //! the boost graph that is encoding the graph grammar rule
		CoreGraph core;

		  //! the id or description of this rules
		std::string id;

		  //! the wildcard string to use for left side pattern matching
		std::string * wildcard;


	public:
		/////////////////  CONSISTENCY TYPEDEFS  /////////////////////////

		  //! consistency code : everything fine
		using Rule_ConCodes::C_Consistent;
		  //! consistency code : no rule ID is present
		using Rule_ConCodes::C_NoID;

		///////////////////////////////////////////////////////////////////////

		  //! A container that stores the information of a context side of 
		  //! a graph grammar rule like left, right, and context.
		class RuleSide {
		protected:
			
			  // to enable direct access from enclosing class
			friend class Rule;

			  /*! Construction
			   * @param nodes the node indices in the core graph that are 
			   *        present in this Rule side
			   * @param constraints the additional constraint for this rule side
			   *        where the node IDs are given according the core of the
			   *        rule
			   */
			RuleSide(	const NodeVec& nodes
						, const MatchConstraintVec * const constraints = NULL );
			
			  //! copy construction
			  //! @param toCopy the object to make a copy of
			RuleSide(	const RuleSide& toCopy );

			  //! destruction
			~RuleSide();

			  //! assignment operator
			  //! @param toCopy the object to make this a copy of
			  //! @return access to the changed *this object
			RuleSide&
			operator=(	const RuleSide& toCopy );

		public:
			
			  //! the list of node indices in the core graph present in this 
			  //! RuleSide
			NodeVec nodes;
			  //! the additional constraints for this RuleSide
			MatchConstraintVec constraints;
		};
		
		///////////////////////////////////////////////////////////////////////


		 /*!
		  * Class to define copy-and-Paste operations for adjacent edges of
		  * nodes to be removed. They are used to relink dangling edges to
		  * other nodes of the result graph.
		  */
		class RuleCnP {

		public:

			  //! the type that defines a set of labels
			typedef std::set< std::string > LabelSet;

			 /*! Default construction initializting with INT_MAX
			  */
			RuleCnP( );

			 /*! Construction
			  * @param source the node that will be deleted and which is the
			  *               source of the dangling edges
			  * @param pasteID the node where the dangling edges will be
			  *               reconnected to
			  */
			RuleCnP( const size_t source
					, const size_t pasteID );

			 /*! Construction
			  * @param source the node that will be deleted and which is the
			  *               source of the dangling edges
			  * @param pasteID the node where the dangling edges will be
			  *               reconnected to
			  * @param target the target node id of the edges to copy-and-paste
			  */
			RuleCnP( const size_t source
					, const size_t pasteID
					, const size_t target );

			 /*! Construction
			  * @param source the node that will be deleted and which is the
			  *               source of the dangling edges
			  * @param pasteID the node where the dangling edges will be
			  *               reconnected to
			  * @param edgeLabels the labels of edges that have to be
			  *   copy-and-pasted; if the set is empty, all edges will be copied.
			  */
			RuleCnP( const size_t source
					, const size_t pasteID
					, const LabelSet & edgeLabels );

			 /*! Construction
			  * @param source the node that will be deleted and which is the
			  *               source of the dangling edges
			  * @param pasteID the node where the dangling edges will be
			  *               reconnected to
			  * @param target the target node id of the edges to copy-and-paste
			  * @param edgeLabels the labels of edges that have to be
			  *   copy-and-pasted; if the set is empty, all edges will be copied.
			  */
			RuleCnP( const size_t source
					, const size_t pasteID
					, const size_t target
					, const LabelSet & edgeLabels );

			  //! the node that will be deleted and which is the source of
			  //! the dangling edges
			size_t source;
			  //! the node where the dangling edges will be reconnected to
			size_t pasteID;
			  //! optional: the target id of the edges to copy-and-paste
			size_t target;
			  //! the labels of edges that have to be copy-and-pasted;
			  //! if the set is empty, all edges will be copied.
			LabelSet edgeLabels;

			  /*!
			   * Order comparison: this object is smaller iff
			   * (this.source < larger.source)
			   * or
			   * (this.source == larger.source and this.target < larger.target).
			   * @param larger the copy-and-paste object to compare to
			   * @return whether this object is of smaller order or not.
			   */
			bool
			operator<( const RuleCnP& larger ) const;

		};

		  //! container definition to store copy-and-Paste operations to do
		  //! the key value defines the source node of the copy-and-Paste
		typedef std::map< size_t, std::set<RuleCnP> > CopyAndPasteOperations;

		///////////////////////////////////////////////////////////////////////

		  //! Construction of a rule based on the given boost graph based core
		  //! description.
		  //! @param core the boost graph that describes the graph grammar rule
		  //! @param id the identifier or description of the graph grammar rule
		  //! @param wildcardToUse the wildcard string to be used for the left
		  //!        side pattern matching, or NULL if no wildcard to be applied
		Rule (	const CoreGraph & core
				, const std::string& id
				, const std::string* wildcardToUse = NULL );
		
		  //! Construction of a rule based on the given boost graph based core
		  //! description and the additional constraints of the left side
		  //! matching.
		  //! @param core the boost graph that describes the graph grammar rule
		  //! @param id the identifier or description of the graph grammar rule
		  //! @param matchConstraints the additional constraints for the
		  //!        matching of the left side of the rule
		  //! @param wildcardToUse the wildcard string to be used for the left
		  //!        side pattern matching, or NULL if no wildcard to be applied
		Rule (	const CoreGraph & core
				, const std::string& id
				, const MatchConstraintVec & matchConstraints
				, const std::string* wildcardToUse = NULL );

		  //! Construction of a rule based on the given boost graph based core
		  //! description and the additional constraints of the left side
		  //! matching.
		  //! @param core the boost graph that describes the graph grammar rule
		  //! @param id the identifier or description of the graph grammar rule
		  //! @param matchConstraints the additional constraints for the
		  //!        matching of the left side of the rule
		  //! @param copyAndPaste the copy-and-Paste operations to perform
		  //! @param wildcardToUse the wildcard string to be used for the left
		  //!        side pattern matching, or NULL if no wildcard to be applied
		Rule (	const CoreGraph & core
				, const std::string& id
				, const MatchConstraintVec & matchConstraints
				, const CopyAndPasteOperations & copyAndPaste
				, const std::string* wildcardToUse = NULL );

		  //! Copy construction.
		  //! @param toCopy the Rule to make this a copy of
		Rule ( const Rule & toCopy );
		
		  //! Assignment operator
		  //! @param toCopy the Rule to make this a copy of
		  //! @return *this
		Rule &
		operator =( const Rule & toCopy );
		
		  //! Destruction
		virtual 
		~Rule();
		
		
		  /*!
		   * Access to the core graph that encodes the whole Rule.
		   * @return the Rule's core graph
		   */
		const CoreGraph&
		getCore(void) const;

	protected:
	
		  //! the information describing the left side pattern of the Rule
		RuleSide LeftSide;
		  //! the information describing the right side (result) of the Rule
		RuleSide RightSide;
		  //! the information describing the invariant context of the Rule
		RuleSide Context;
		
		  //! the copy-and-Paste operations to perform
		CopyAndPasteOperations copyAndPaste;

		  /*! Utility member function to derive the node information 
		   * that describes a certain part of the Rule.
		   * @param graph the core graph that contains the whole information
		   * @param level identifies what RuleSide is of interest
		   * @return the vector of node indices of the requested RuleSide
		   */
		static 
		NodeVec 
		getRuleSideNodes(	const CoreGraph &graph,
							const RuleContext &level );
		
	public:


		  /*! Checks whether or not this rule is a consistent one. If not, a
		   * combination of according error codes is given.
		   * @return the value is either C_Consistent, or a product of all other
		   *         consistency error codes encountered
		   */
		virtual
		size_t
		isConsistent( void ) const;
		
		  /*!
		   * Writes a description of the consistency status or errors, encoded
		   * in a consistency code produced by isConsistent*(...), to a given
		   * outstream. The function returns whether or not an error occured.
		   *
		   * @param consistencyCode the error code to parse, produced by
		   *           a call to isConsistent(...)
		   * @param errorStream the output stream to write the error decription
		   *           to
		   * @param completeCheck if true: tries to decode the whole
		   *           consistencyCode and reports an error message if this is
		   *           not possible. if false: decodes and reports only the
		   *           known error codes from consistencyCode
		   * @return true if no error is encoded; false otherwise
		   */
		virtual
		bool
		decodeConsistencyStatus(	const size_t consistencyCode
									, std::ostream& errorStream
									, const bool completeCheck = true ) ;

		  /*!
		   * Access to the left side pattern of the rule that has to be matched
		   * when the Rule is applied.
		   * @return the left side pattern
		   */
		const RuleSide&
		getLeftSide(void) const;

		  /*!
		   * Access to the right side pattern that is produced after the Rule
		   * was applied.
		   * @return the right side pattern
		   */
		const RuleSide&
		getRightSide(void) const;

		  /*!
		   * Access to the context of the Rule that is fixed during the Rule
		   * application and that is present in the left and right side of the
		   * Rule.
		   * @return the context
		   */
		const RuleSide&
		getContext(void) const;
		
		  /*!
		   * Access to the string identifier or description of the Rule.
		   * @return the identifier
		   */ 
		const std::string&
		getID(void) const;
		
		  /*!
		   * Sets the string identifier or description of the Rule.
		   * @param newID the new identifier to set
		   */
		void
		setID( const std::string& newID );
		
		  /*!
		   * Access to the wildcard string to use for left side pattern matching
		   * @return the wildcard or NULL if no wildcard is to be applied
		   */
		virtual
		const std::string*
		getUsedWildcard(void) const;

		  /*!
		   * Sets the wildcard string to use for left side pattern matching.
		   * To disable wildcard matching set the wildcard to NULL.
		   * @param wildcardToUse the wildcard to be applied or NULL if no
		   *        wildcard matching should be done.
		   */
		virtual
		void
		setUsedWildcard( const std::string* const wildcardToUse );

		  /*!
		   * Access to the match constraints to fulfill for each left side match
		   * @return the match constraint container
		   */
		const MatchConstraintVec &
		getMatchConstraints( void ) const;

		  /*!
		   * Access to the copy-and-Paste operations to perform
		   * @return the CnP operation container
		   */
		const CopyAndPasteOperations &
		getCopyAndPasteOperations( void ) const;


	protected:

		  /*! Checks the context assumptions of all edges, ie.
		   * - context edges have to connect context nodes,
		   * - left side edges have to connect left or context nodes, and
		   * - right side edges have to connect right or context nodes.
		   *
		   * @param coreGraph the rule core graph to check
		   * @return true if all edges show correct context use; false otherwise
		   */
		static
		bool
		checkEdgeContext( const CoreGraph & coreGraph );


	};
	
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}

  // implementation
#include "ggl/Rule.icc"

#endif /*RULE_HH_*/
