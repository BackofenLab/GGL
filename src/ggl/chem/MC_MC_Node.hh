#ifndef GGL_CHEM_MC_MC_NODE_HH_
#define GGL_CHEM_MC_MC_NODE_HH_

#include <sgm/MC_Node.hh>

#include <ggl/chem/MoleculeUtil.hh>


namespace ggl {
namespace chem {



//=========================================================================

	  /*! @brief Match constraints of MoleculeComponents
	   *
	   * This an sgm::MC_NodeAdjacency constraint for MoleculeComponent
	   * descriptions. Therefore, it assumes valid (complex) atom labels for
	   * all nodes.
	   *
	   * NOTE: Only the atom identifier is used for comparisons. All charge,
	   * class, or proton information is ignored.
	   *
	   * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */
	class MC_MC_NodeAdjacency : public sgm::MC_NodeAdjacency {


		  //! the type that defines a set of labels
		typedef MC_NodeAdjacency::LabelSet LabelSet;

	public:

		  //! the index of the node this constraint is for
		using MC_NodeAdjacency::constrainedNodeID;

		  //! the relational operator to be applied using 'count'
		using MC_NodeAdjacency::op;

		  //! the number of adjacent nodes/edges fulfilling the label condition
		using MC_NodeAdjacency::count;

		  //! the node labels of adjacent nodes that have to be counted
		using MC_NodeAdjacency::nodeLabels;

		  //! the labels of adjacent edges that have to be counted
		using MC_NodeAdjacency::edgeLabels;

		  /*!
		   * Default construction of the match constraint : undefined values are
		   * initialized with INT_MAX
		   */
		MC_MC_NodeAdjacency()
		 : sgm::MC_NodeAdjacency()
		{}

		  /*!
		   * Copy construction of the match constraint
		   */
		MC_MC_NodeAdjacency( const MC_NodeAdjacency& toCopy )
		 : MC_NodeAdjacency(toCopy)
		{}

		  /*!
		   * Copy construction of the match constraint
		   */
		MC_MC_NodeAdjacency( const MC_MC_NodeAdjacency& toCopy )
		 : sgm::MC_NodeAdjacency(toCopy)
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
		MC_MC_NodeAdjacency(	const size_t constrainedNodeID
								, const MC_Operator op
								, const size_t count
								, const LabelSet & nodeLabels
								, const LabelSet & edgeLabels
						)
		 :  sgm::MC_NodeAdjacency( constrainedNodeID , op , count , nodeLabels , edgeLabels )
		{}

		  /*!
		   * Construction of the match constraint
		   *
		   * @param constrainedNodeID the index of the node this constraint is for
		   * @param op the relational operator to be applied using 'count'
		   * @param count the number of adjacent nodes/edges fulfilling the label condition
		   * @param nodeLabel a single node label of adjacent nodes that has to be counted
		   */
		MC_MC_NodeAdjacency(	const size_t constrainedNodeID
								, const MC_Operator op
								, const size_t count
								, const std::string & nodeLabel )
		 : sgm::MC_NodeAdjacency( constrainedNodeID, op, count, nodeLabel )
		{}

		  //! destruction
		virtual
		~MC_MC_NodeAdjacency()
		{}


		 /*!
		  * Checks whether or not a match on a given target fulfills the
		  * additional node adjacency constraint for the pattern matching.
		  *
		  * NOTE: If node labels are given, only checking the atom identifier,
		  * not the whole node label is checked against the label set.
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
		isValidMatch(	const sgm::Pattern_Interface & pattern,
						const sgm::Graph_Interface & target,
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
					std::unordered_map<sgm::Graph_Interface::IndexType, char>
#elif HAVE_TR1_UNORDERED_MAP > 0
					std::tr1::unordered_map<sgm::Graph_Interface::IndexType, char>
#elif HAVE_GNU_HASH_MAP > 0
					__gnu_cxx::hash_map<sgm::Graph_Interface::IndexType, char>
#else
					std::map<sgm::Graph_Interface::IndexType, char>
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
						if (nodeWildcard || nodeLabels.find(MoleculeUtil::getAtom(target.getNodeLabel(curEdge->getToIndex()))) != nodeLabels.end()) {
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
							&& (allNodes || nodeLabels.find(MoleculeUtil::getAtom(target.getNodeLabel(curEdge->getToIndex()))) != nodeLabels.end()))
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
//
//			  // will hold the number of label matches among adjacent nodes in target graph
//			size_t hits = 0;
//
//			  // check if nodeLabels && edgeLabels empty --> placeholder for matching anything
//			if (nodeLabels.size() == 0 && edgeLabels.size() == 0 ) {
//				  // count adjacent edges
//				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
//						curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
//					curEdge != curEdgeEnd; ++curEdge )
//				{
//					hits++;
//				}
//			} else
//				// ONLY NODE LABELS SPECIFIED
//			if (edgeLabels.size() == 0) {
//				  // iterate through all adjacent nodes of match[rC.constrainedNodeID]
//				  // and check node label
//				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
//						curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
//					curEdge != curEdgeEnd; ++curEdge )
//				{
//					  // check if atom label of the end of the current edge is among
//					  // the labels to count
//					if (nodeLabels.find(MoleculeUtil::getAtom(target.getNodeLabel(curEdge->getToIndex()))) != nodeLabels.end()) {
//						hits++;
//					}
//				}
//			} else
//				// ONLY EDGE LABELS SPECIFIED
//			if (nodeLabels.size() == 0) {
//				  // iterate through all edges of match[constrainedNodeID]
//				  // and check edge label
//				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
//						curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
//					curEdge != curEdgeEnd; ++curEdge )
//				{
//					  // check if edge label is among the labels to count
//					if (edgeLabels.find(curEdge->getEdgeLabel()) != edgeLabels.end()) {
//						hits++;
//					}
//				}
//				// BOTH : NODE AND EDGE LABELS SPECIFIED
//			} else {
//				  // iterate through all edges of match[constrainedNodeID]
//				  // and check edge AND node label
//				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(matchedTargetID),
//						curEdgeEnd = target.getOutEdgesEnd(matchedTargetID);
//					curEdge != curEdgeEnd; ++curEdge )
//				{
//					  // check if edge and node label are among the labels to count
//					if (edgeLabels.find(curEdge->getEdgeLabel()) != edgeLabels.end()
//						&& nodeLabels.find(MoleculeUtil::getAtom(target.getNodeLabel(curEdge->getToIndex()))) != nodeLabels.end())
//					{
//						hits++;
//					}
//				}
//			}
//
//			  // check if constraint is fulfilled based on relational operator
//			switch( op ) {
//				case MC_EQ : return ( hits == count ); break;
//				case MC_L  : return ( hits  < count ); break;
//				case MC_G  : return ( hits  > count ); break;
//				case MC_LQ : return ( hits <= count ); break;
//				case MC_GQ : return ( hits >= count ); break;
//				default : assert(false /* UNSUPPORTED RELATIONAL OPERATOR */);
//			}
//
//			  // final decision, should never be called
//			return false;
		}

		 /*!
		  * Creates a new MC_NodeAdjacency heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_NodeAdjacency object that equals *this
		  */
		virtual
		MC_MC_NodeAdjacency*
		clone( void ) const
		{
			return new MC_MC_NodeAdjacency(*this);
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
		MC_MC_NodeAdjacency*
		remap( const sgm::Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
		{
			assert(this->constrainedNodeID < old2newIndexMapping.size());
			  // check if this node is an unmatched node and thus to be ignored
			if (old2newIndexMapping.at(this->constrainedNodeID)==unmatchedIndex) {
				return NULL;
			}
			  // create copy
			MC_MC_NodeAdjacency* copy = new MC_MC_NodeAdjacency(*this);
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
		operator==( const MC_MC_NodeAdjacency& toCompare ) const
		{
			return sgm::MC_NodeAdjacency::operator==(toCompare);
		}


	};

//=========================================================================


	 /*!
	  * A an sgm::MC_NodeLabel constraint for MoleculeComponent
	  * descriptions. Therefore, it assumes valid (complex) atom labels for
	  * all nodes.
	  *
	  * NOTE: Only the atom identifier is used for comparisons. All charge,
	  * class, or proton information is ignored.
	  *
	  * NOTE: This constraint is only useful if the according pattern node shows
	  *       a wildcard label such that it can be matched on any target node.
	  *
	  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MC_MC_NodeLabel : public sgm::MC_NodeLabel {

	public:

		  //! the type that defines a set of labels
		typedef MC_NodeLabel::LabelSet LabelSet;


	public:

		  //! the node ID to be constrained
		using MC_NodeLabel::constrainedNodeID;

		  //! the set of labels this node can be mapped on
		  //! NOTE: if empty, no matching will be allowed !!!
		using MC_NodeLabel::nodeLabels;

		  //! the type of comparison to be applied, i.e. if to ensure that the
		  //! edge label is among the edgeLabels (ALLOWED) or that it is not
		  //! present (FORBIDDEN)
		using MC_NodeLabel::compareType;

		  /*!
		   * Default construction of the match constraint : undefined values are
		   * initialized with INT_MAX or empty sets
		   */
		MC_MC_NodeLabel()
		 : MC_NodeLabel()
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedNodeID the node ID to be constrained
		MC_MC_NodeLabel( const size_t constrainedNodeID )
		  :	MC_NodeLabel(constrainedNodeID)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedNodeID the node ID to be constrained
		  //! @param nodeLabels the allowed node labels to be matched on
		MC_MC_NodeLabel( const size_t constrainedNodeID
						, const LabelSet& nodeLabels)
		  :	MC_NodeLabel( constrainedNodeID, nodeLabels )
		{}


		  //! construction where the allowed node labels are empty
		  //! @param constrainedNodeID the node ID to be constrained
		  //! @param nodeLabels the allowed node labels to be matched on
		  //! @param compareType the type of comparison to be applied
		MC_MC_NodeLabel( const size_t constrainedNodeID
						, const LabelSet& nodeLabels
						, const CompareType& compareType)
		  :	MC_NodeLabel( constrainedNodeID, nodeLabels, compareType )
		{}

		  //! copy construction
		  //! @param toCopy the MC_NodeLabel object to copy
		MC_MC_NodeLabel( const MC_NodeLabel& toCopy )
		  :	MC_NodeLabel(toCopy)
		{}

		  //! copy construction
		  //! @param toCopy the MC_MC_NodeLabel object to copy
		MC_MC_NodeLabel( const MC_MC_NodeLabel& toCopy )
		  :	MC_NodeLabel(toCopy)
		{}

		  //! destruction
		virtual
		~MC_MC_NodeLabel()
		{}

		 /*!
		  * Checks whether or not the matched node holds one of the allowed
		  * labels.
		  *
		  * NOTE: If node labels are given, only checking the atom identifier,
		  * not the whole node label is checked against the label set.
		  *
		  * @param pattern the pattern graph that was matched
		  * @param target the target graph the pattern was matched on
		  * @param matchedTargetID matched node index within the target graph
		  * @return true if the match is valid; false if the constraint is
		  *         violated
		  */
		virtual
		bool
		isValidMatch(	const sgm::Pattern_Interface & pattern,
						const sgm::Graph_Interface & target,
						const size_t matchedTargetID ) const
		{
			  // check if the mapped node shows an allowed label
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
				return nodeLabels.find(MoleculeUtil::getAtom(target.getNodeLabel(matchedTargetID)))
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
				return nodeLabels.find(MoleculeUtil::getAtom(target.getNodeLabel(matchedTargetID)))
						== nodeLabels.end();

			default :
				assert(false); /*should never happen*/
			}

			  // if not satisfied so far, wont be satisfied anymore.. ;)
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
		MC_MC_NodeLabel*
		clone( void ) const
		{
			return new MC_MC_NodeLabel(*this);
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
		MC_MC_NodeLabel*
		remap( const sgm::Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
		{
			assert(this->constrainedNodeID < old2newIndexMapping.size());
			  // check if this node is an unmatched node and thus to be ignored
			if (old2newIndexMapping.at(this->constrainedNodeID)==unmatchedIndex) {
				return NULL;
			}
			  // create copy
			MC_MC_NodeLabel* copy = new MC_MC_NodeLabel(*this);
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
		operator==( const MC_MC_NodeLabel& toCompare ) const
		{
			return MC_NodeLabel::operator==(toCompare);
		}



	};

//=========================================================================

}} // namespace ggl


#endif /* GGL_CHEM_MC_MC_NODE_HH_ */

