#ifndef SGM_MC_NODE_HH_
#define SGM_MC_NODE_HH_

#include <string>
#include <set>
#include <climits>
#include <cassert>

#include <sgm/Pattern.hh>
#include <sgm/Graph_Interface.hh>
#include <sgm/Match_Reporter.hh>


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


namespace sgm {


//=========================================================================


	 /*! @brief Interface node match constraints
	  *
	  * A match node constraint describes additional properties that have to be
	  * fulfilled by a given matched node on a given target.
	  *
	  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MC_Node : public Pattern_Interface::Match_Constraint {

	public:

		  //! the node ID to be constrained
		size_t constrainedNodeID;

		  //! construction
		  //! @param constrainedNodeID the node ID to be constrained
		MC_Node( const size_t constrainedNodeID )
		  :	constrainedNodeID(constrainedNodeID)
		{}

		  //! copy construction
		  //! @param toCopy the MC_Node object to copy
		MC_Node( const MC_Node& toCopy )
		  :	constrainedNodeID(toCopy.constrainedNodeID)
		{}

		  //! destruction
		virtual
		~MC_Node()
		{}

		 /*!
		  * Checks whether or not this constraint covers the node with the
		  * given ID.
		  *
		  * @param nodeID the ID of the node of interest
		  * @return true if the node is covered by the constraint; false
		  *         otherwise
		  */
		virtual
		bool
		isConstraining(	const size_t nodeID ) const
		{
			  // check if this ID is constrained
			return nodeID == constrainedNodeID;
		}

		 /*!
		  * Creates a new MC_Node heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_Node object that equals this
		  */
		virtual
		MC_Node*
		clone( void ) const = 0;

		 /*!
		  * Creates a new MC_Node heap object that equals the current
		  * object but uses the new indices given by old2newIndexMapping.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @param old2newIndexMapping the index mapping to be used for the
		  *        remapping
		  * @param unmatchedIndex an optional specific index that marks
		  *        unmatched nodes within old2newIndexMapping. if this
		  *        constrains one of these nodes, no remapping is done and
		  *        NULL is returned
		  * @return a new allocated MC_Node object
		  */
		virtual
		MC_Node*
		remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX ) = 0;


		 /*!
		  * Checks whether or not a match on a given target fulfills the
		  * additional node constraint for the pattern matching.
		  *
		  * @param pattern the pattern graph that was matched
		  * @param target the target graph the pattern was matched on
		  * @param matchedTargetID the matched node index within
		  *        the target graph
		  * @return true if the match is valid; false if the constraint is
		  *         violated
		  */
		virtual
		bool
		isValidMatch(	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						const size_t matchedTargetID ) const = 0;


		 /*!
		  * Checks whether or not a match on a given target fulfills the
		  * additional node constraint for the pattern matching.
		  *
		  * @param pattern the pattern graph that was matched
		  * @param target the target graph the pattern was matched on
		  * @param match the match information for the left side pattern of the
		  *        pattern on the target graph
		  * @return true if the match is valid; false if the constraint is
		  *         violated
		  */
		virtual
		bool
		isValidMatch(	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						const Match & match ) const
		{
			  // forward call to specific function
			return isValidMatch( pattern, target, match.at(constrainedNodeID));
		}
		 /*!
		  * Equality comparison to another match constraint.
		  * @param toCompare the constraint to compare to
		  * @return true if the constraints are equal; false otherwise
		  */
		virtual
		bool
		operator==( const sgm::Pattern_Interface::Match_Constraint & toCompare ) const {
			const MC_Node* pc = dynamic_cast< const MC_Node* >(&toCompare);
			if (pc != NULL) {
				return this->constrainedNodeID == pc->constrainedNodeID;
			} else {
				return false;
			}
		}


	};

//=========================================================================

	  /*! @brief Constrains a node's adjacency for a match
	   *
	   * This match constraint defines the relation of the number of adjacent
	   * nodes/edges that have a certain label (combination).
	   *
	   * The constraint is fulfilled if:
	   *
	   * (N(constrainedNodeID,nodeLabels,edgeLabels) op count)
	   *
	   * where N(constrainedNodeID,nodeLabels,edgeLabels) gives the number of adjacent
	   * edges of constrainedNodeID that show a label within edgeLabels and where the
	   * connected node shows on of nodeLabels.
	   *
	   * current idea of GML encoding:
	   *
	   * \verbatim
*****************************************************
 constrainAdj [
   id   ::= INTEGER
   op     ::= '=' | '<' | '>'
   count  ::= DIGIT+
   nodeLabels [
     label ::= STRING
   ]
   edgeLabels [
     label ::= STRING
   ]
 ]
*****************************************************
	   * \endverbatim
	   *
	   * a missing nodeLabels/edgeLabels is used to match anything.
	   *
	   * Examples:
	   *
	   * \verbatim
*****************************************************
 constrainAdj [
   id 5
   op <
   count 2
   nodeLabels [ "N1" "N2" ]
   edgeLabels [ "E" ]
 ]
*****************************************************
	   * \endverbatim
	   *
	   * constrains the 5th node to show less than 2 edges with label "E" that
	   * target a node with label "N1" or "N2".
	   *
	   * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */
	class MC_NodeAdjacency : public MC_Node {

	public:

		  /*!
		   * Possible operators to be applied onto the count of matched adjacent
		   * edges and nodes.
		   */
		enum MC_Operator {	MC_EQ	//!< equality operator
							, MC_L	//!< less operator
							, MC_G	//!< greater operator
							, MC_LQ	//!< less or equal operator
							, MC_GQ	//!< greater or equal operator
		};

		  //! the type that defines a set of labels
		typedef std::set< std::string > LabelSet;

	public:

		  //! the index of the node this constraint is for
		using MC_Node::constrainedNodeID;

		  //! the relational operator to be applied using 'count'
		MC_Operator op;

		  //! the number of adjacent nodes/edges fulfilling the label condition
		size_t count;

		  //! the node labels of adjacent nodes that have to be counted
		LabelSet nodeLabels;

		  //! the labels of adjacent edges that have to be counted
		LabelSet edgeLabels;

		  /*!
		   * Default construction of the match constraint : undefined values are
		   * initialized with INT_MAX
		   */
		MC_NodeAdjacency()
		 : MC_Node((size_t)INT_MAX)
		   , op(MC_EQ)
		   , count((size_t)INT_MAX)
		   , nodeLabels()
		   , edgeLabels()
		{}

		  /*!
		   * Copy construction of the match constraint
		   * @param toCopy the object to make this a copy of
		   */
		MC_NodeAdjacency( const MC_NodeAdjacency& toCopy )
		 : MC_Node( toCopy )
		   , op(toCopy.op)
		   , count(toCopy.count)
		   , nodeLabels(toCopy.nodeLabels)
		   , edgeLabels(toCopy.edgeLabels)
		{}

		  /*!
		   * Construction of the match constraint
		   *
		   * @param constrainedNodeID the index of the node this constraint is for
		   * @param op the relational operator to be applied using 'count'
		   * @param count the number of adjacent nodes/edges fulfilling the label condition
		   * @param nodeLabels the node labels of adjacent nodes that have to be counted
		   * @param edgeLabels the labels of adjacent edges that have to be counted
		   */
		MC_NodeAdjacency(	const size_t constrainedNodeID
						, const MC_Operator op
						, const size_t count
						, const LabelSet & nodeLabels
						, const LabelSet & edgeLabels
						)
		 : MC_Node(constrainedNodeID)
		   , op(op)
		   , count(count)
		   , nodeLabels(nodeLabels)
		   , edgeLabels(edgeLabels)
		{}

		  /*!
		   * Construction of the match constraint
		   *
		   * @param constrainedNodeID the index of the node this constraint is for
		   * @param op the relational operator to be applied using 'count'
		   * @param count the number of adjacent nodes/edges fulfilling the label condition
		   * @param nodeLabel a single node label of adjacent nodes that has to be counted
		   */
		MC_NodeAdjacency(	const size_t constrainedNodeID
						, const MC_Operator op
						, const size_t count
						, const std::string & nodeLabel )
		 : MC_Node(constrainedNodeID)
		   , op(op)
		   , count(count)
		   , nodeLabels()
		   , edgeLabels()
		{
			nodeLabels.insert(nodeLabel);
		}

		  //! destruction
		virtual
		~MC_NodeAdjacency()
		{}

		 /*!
		  * Checks whether or not a match on a given target fulfills the
		  * additional node adjacency constraint for the pattern matching.
		  *
		  * @param pattern the pattern graph that was matched
		  * @param target the target graph the pattern was matched on
		  * @param matchedTargetID the matched node index within
		  *        the target graph
		  * @return true if the match is valid; false if the constraint is
		  *         violated
		  */
		virtual
		bool
		isValidMatch(	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						const size_t matchedTargetID ) const
		{
			  // will hold the number of label matches among adjacent nodes in target graph
			size_t hits = 0;

			  // check if wildcard to use and among allowed node labels
			const bool nodeWildcard = pattern.getUsedWildcard() != NULL
									&& nodeLabels.find(*(pattern.getUsedWildcard())) != nodeLabels.end();
			  // check if wildcard to use and among allowed edge labels
			const bool edgeWildcard = pattern.getUsedWildcard() != NULL
									&& edgeLabels.find(*(pattern.getUsedWildcard())) != edgeLabels.end();

			  // special case of checking only adjacent nodes, ignoring the edges
			  // -> no edge but only node labels specified but
			if (edgeLabels.empty() && !(nodeLabels.empty())) {
				// ONLY NODE LABELS SPECIFIED
				  // will hold the observed adjacent indices to enable correct
				  // node label counting in case multiple parallel edges are
				  // present between this and an adjacent node
#if HAVE_UNORDERED_MAP > 0
					std::unordered_map<Graph_Interface::IndexType, char>
#elif HAVE_TR1_UNORDERED_MAP > 0
					std::tr1::unordered_map<Graph_Interface::IndexType, char>
#elif HAVE_GNU_HASH_MAP > 0
					__gnu_cxx::hash_map<Graph_Interface::IndexType, char>
#else
					std::map<Graph_Interface::IndexType, char>
#endif
						 observedNodes;
				  // iterate through all adjacent nodes of match[rC.constrainedNodeID]
				  // and check node label
				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
						curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
					curEdge != curEdgeEnd; ++curEdge )
				{
					  // check if the current target node was already checked
					if (observedNodes.find(curEdge->getToIndex()) == observedNodes.end()) {
						  // add to observed indices
						observedNodes[curEdge->getToIndex()] = 'T';
						  // check if label of the end of the current edge is among
						  // the labels to count
						if (nodeWildcard || nodeLabels.find(target.getNodeLabel(curEdge->getToIndex())) != nodeLabels.end()) {
							hits++;
						}
					}
				}
			} else {

				  // check if all adjacent nodes are to be counted
				const bool allNodes = nodeLabels.empty() || nodeWildcard;
				  // check if all edges are to be counted
				const bool allEdges = edgeLabels.empty() || edgeWildcard;

				  // check if all edges to be counted --> placeholder for matching anything
				if ( allNodes && allEdges )
				{
					  // count adjacent edges
					for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
							curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
						curEdge != curEdgeEnd; ++curEdge )
					{
						hits++;
					}
				} else {
					  // iterate through all edges of match[constrainedNodeID]
					  // and check edge AND node label
					for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
							curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
						curEdge != curEdgeEnd; ++curEdge )
					{
						  // check if edge and node label are among the labels to count
						if ((allEdges || edgeLabels.find(curEdge->getEdgeLabel()) != edgeLabels.end())
							&& (allNodes || nodeLabels.find(target.getNodeLabel(curEdge->getToIndex())) != nodeLabels.end()))
						{
							hits++;
						}
					}
				}
			}
			  // check if constraint is fulfilled based on relational operator
			switch( op ) {
				case MC_EQ : return ( hits == count ); break;
				case MC_L  : return ( hits  < count ); break;
				case MC_G  : return ( hits  > count ); break;
				case MC_LQ : return ( hits <= count ); break;
				case MC_GQ : return ( hits >= count ); break;
				default : assert(false /* UNSUPPORTED RELATIONAL OPERATOR */);
			}

			  // final decision, should never be called
			return false;
		}

		 /*!
		  * Creates a new MC_NodeAdjacency heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_NodeAdjacency object that equals *this
		  */
		virtual
		MC_NodeAdjacency*
		clone( void ) const
		{
			return new MC_NodeAdjacency(*this);
		}

		virtual
		bool
		isConstrainedLabel(	const std::string & label ) const
		{
			  // check if among the node labels
			for (LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end(); ++l)
				if (l->compare(label)==0)
					return true;
			  // check if among the edge labels
			for (LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end(); ++l)
				if (l->compare(label)==0)
					return true;
			  // not found
			return false;
		}


		 /*!
		  * Creates a new MC_NodeAdjacency heap object that equals the current
		  * object but uses the new indices given by old2newIndexMapping.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @param old2newIndexMapping the index mapping to be used for the
		  *        remapping
		  * @param unmatchedIndex an optional specific index that marks
		  *        unmatched nodes within old2newIndexMapping. if this
		  *        constrains one of these nodes, no remapping is done and
		  *        NULL is returned
		  * @return a new allocated MC_NodeAdjacency object
		  */
		virtual
		MC_NodeAdjacency*
		remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
		{
			assert(this->constrainedNodeID < old2newIndexMapping.size());

			  // check if this constrains an unmatched node
			if (old2newIndexMapping.at(this->constrainedNodeID) == unmatchedIndex)
			{
				return NULL;
			}

			  // create copy
			MC_NodeAdjacency* copy = new MC_NodeAdjacency(*this);
			  // do remapping
			copy->constrainedNodeID = old2newIndexMapping.at(copy->constrainedNodeID);
			  // return remapped copy
			return copy;
		}


		 /*!
		  * Equality comparison to another match constraint.
		  * @param toCompare the constraint to compare to
		  * @return true if the constraints are equal; false otherwise
		  */
		virtual
		bool
		operator==( const MC_NodeAdjacency& toCompare ) const
		{
			const MC_NodeAdjacency* pc = dynamic_cast< const MC_NodeAdjacency* >(&toCompare);
			if (pc != NULL) {
				bool isEqual =
					MC_Node::operator==( toCompare )
					&& this->op == pc->op
					&& this->count == pc->count
					&& this->nodeLabels.size() == pc->nodeLabels.size()
					&& this->edgeLabels.size() == pc->edgeLabels.size()
					;
				  // check if all node labels present
				for (LabelSet::const_iterator l=this->nodeLabels.begin();
						isEqual && l!=this->nodeLabels.end(); ++l)
				{
					isEqual = pc->nodeLabels.find(*l)!=pc->nodeLabels.end();
				}
				  // check if all edge labels present
				for (LabelSet::const_iterator l=this->edgeLabels.begin();
						isEqual && l!=this->edgeLabels.end(); ++l)
				{
					isEqual = pc->edgeLabels.find(*l)!=pc->edgeLabels.end();
				}
				  // return final comparison result
				return isEqual;
			} else {
				return false;
			}
		}


	};

//=========================================================================


	 /*! @brief Constrains a node's labels for a matching
	  *
	  * A match constraint that restricts the allowed matched labels for a given
	  * matched node to a set of maintained labels. Depending on the operator
	  * type, the maintained node label set describes the set of allowed labels
	  * (op =) or the set of forbidden labels (op !). The operator defaults to
	  * allowed (op =).
	  *
	  * NOTE: This constraint is only useful if the according pattern node shows
	  *       a wildcard label such that it can be matched on any target node.
	  *
	   * current idea of GML encoding:
	   *
	   * \verbatim
*****************************************************
 constrainNode [
   id   ::= INTEGER
   op   ::= OPERATOR { = | ! }
   nodeLabels [
     label ::= STRING
   ]
 ]
*****************************************************
	   * \endverbatim
	   *
	   * a missing nodeLabels/edgeLabels is used to match anything.
	   *
	   * Examples:
	   *
	   * \verbatim
*****************************************************
 constrainNode [
   id 1
   nodeLabels [
     label "C"
     label "P"
   ]
 ]
*****************************************************
	   * \endverbatim
	   *
	   * constrains the allowed labels of the node with index 1 to "C" or "P".
	   *
	  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MC_NodeLabel : public MC_Node {

	public:

		  //! the type that defines a set of labels
		typedef std::set< std::string > LabelSet;

		enum CompareType { ALLOWED, FORBIDDEN };


	public:

		  //! the node ID to be constrained
		using MC_Node::constrainedNodeID;

		  //! the set of labels this node can be mapped on
		  //! NOTE: if empty, no matching will be allowed !!!
		LabelSet nodeLabels;

		  //! the type of comparison to be applied, i.e. if to ensure that the
		  //! edge label is among the edgeLabels (ALLOWED) or that it is not
		  //! present (FORBIDDEN)
		CompareType compareType;

		  /*!
		   * Default construction of the match constraint : undefined values are
		   * initialized with INT_MAX or empty sets
		   */
		MC_NodeLabel()
		 : MC_Node((size_t)INT_MAX)
		   , nodeLabels()
			, compareType(ALLOWED)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedNodeID the node ID to be constrained
		MC_NodeLabel( const size_t constrainedNodeID )
		  :	MC_Node(constrainedNodeID)
			, nodeLabels()
			, compareType(ALLOWED)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedNodeID the node ID to be constrained
		  //! @param nodeLabels the allowed node labels to be matched on
		MC_NodeLabel( const size_t constrainedNodeID
						, const LabelSet& nodeLabels)
		  :	MC_Node( constrainedNodeID )
			, nodeLabels( nodeLabels )
			, compareType(ALLOWED)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedNodeID the node ID to be constrained
		  //! @param nodeLabels the allowed node labels to be matched on
		  //! @param compareType the type of comparison to be applied
		MC_NodeLabel( const size_t constrainedNodeID
						, const LabelSet& nodeLabels
						, const CompareType& compareType)
		  :	MC_Node( constrainedNodeID )
			, nodeLabels( nodeLabels )
			, compareType( compareType )
		{}

		  //! copy construction
		  //! @param toCopy the MC_NodeLabel object to copy
		MC_NodeLabel( const MC_NodeLabel& toCopy )
		  :	MC_Node(toCopy)
			, nodeLabels( toCopy.nodeLabels )
			, compareType( toCopy.compareType )
		{}

		  //! destruction
		virtual
		~MC_NodeLabel()
		{}

		 /*!
		  * Checks whether or not the matched node holds one of the allowed
		  * labels.
		  *
		  * @param pattern the pattern graph that was matched
		  * @param target the target graph the pattern was matched on
		  * @param matchedTargetID the matched node index within
		  *        the target graph
		  * @return true if the match is valid; false if the constraint is
		  *         violated
		  */
		virtual
		bool
		isValidMatch(	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						const size_t matchedTargetID ) const
		{
			  // handle comparison according to given comparison type
			switch (compareType) {

			case ALLOWED : // check if all edge labels are found
				  // no node label given -> no allowed label available
				if (nodeLabels.empty()) {
					return false;
				}
				  // check if a wildcard is to be used and is among the allowed labels
				if (pattern.getUsedWildcard() != NULL
					&& nodeLabels.find(*(pattern.getUsedWildcard())) != nodeLabels.end() )
				{
					return true;
				}
				  // check if the mapped node shows an allowed label
				return nodeLabels.find(target.getNodeLabel(matchedTargetID))
						!= nodeLabels.end();

			case FORBIDDEN : // check that all edge labels are NOT found
				  // no node label given -> no forbidden label available
				if (nodeLabels.empty()) {
					return true;
				}
				  // check if a wildcard is to be used and is among the allowed labels
				if (pattern.getUsedWildcard() != NULL
					&& nodeLabels.find(*(pattern.getUsedWildcard())) != nodeLabels.end() )
				{
					return false;
				}
				  // check if the mapped node shows no forbidden label
				return nodeLabels.find(target.getNodeLabel(matchedTargetID))
						== nodeLabels.end();

			default :
				assert(false); /*should never happen*/
			}

			  // if not satisfied so far, wont be satisfied anymore.. ;)
			return false;
		}

		virtual
		bool
		isConstrainedLabel(	const std::string & label ) const
		{
			  // check if among the node labels
			for (LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end(); ++l)
				if (l->compare(label)==0)
					return true;
			  // not found
			return false;
		}


		 /*!
		  * Creates a new MC_NodeLabel heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_NodeLabel object that equals *this
		  */
		virtual
		MC_NodeLabel*
		clone( void ) const
		{
			return new MC_NodeLabel(*this);
		}


		 /*!
		  * Creates a new MC_NodeLabel heap object that equals the current
		  * object but uses the new indices given by old2newIndexMapping.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @param old2newIndexMapping the index mapping to be used for the
		  *        remapping
		  * @param unmatchedIndex an optional specific index that marks
		  *        unmatched nodes within old2newIndexMapping. if this
		  *        constrains one of these nodes, no remapping is done and
		  *        NULL is returned
		  * @return a new allocated MC_NodeLabel object
		  */
		virtual
		MC_NodeLabel*
		remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
		{
			assert(this->constrainedNodeID < old2newIndexMapping.size());

			  // check if this constrains an unmatched node
			if (old2newIndexMapping.at(this->constrainedNodeID) == unmatchedIndex)
			{
				return NULL;
			}

			  // create copy
			MC_NodeLabel* copy = new MC_NodeLabel(*this);
			  // do remapping
			copy->constrainedNodeID = old2newIndexMapping.at(copy->constrainedNodeID);
			  // return remapped copy
			return copy;
		}


		 /*!
		  * Equality comparison to another match constraint.
		  * @param toCompare the constraint to compare to
		  * @return true if the constraints are equal; false otherwise
		  */
		virtual
		bool
		operator==( const MC_NodeLabel& toCompare ) const
		{
			const MC_NodeLabel* pc = dynamic_cast< const MC_NodeLabel* >(&toCompare);
			if (pc != NULL) {
				bool isEqual =
					MC_Node::operator==( toCompare )
					&& this->nodeLabels.size() == pc->nodeLabels.size()
					&& this->compareType == pc->compareType
					;
				  // check if all node labels present
				for (LabelSet::const_iterator l=this->nodeLabels.begin();
						isEqual && l!=this->nodeLabels.end(); ++l)
				{
					isEqual = pc->nodeLabels.find(*l)!=pc->nodeLabels.end();
				}
				  // return final comparison result
				return isEqual;
			} else {
				return false;
			}
		}



	};

//=========================================================================

} // namespace sgm


#endif /* SGM_MC_NODE_HH_ */

