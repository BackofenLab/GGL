#ifndef SGM_MC_EDGE_HH_
#define SGM_MC_EDGE_HH_

#include <string>
#include <set>
#include <climits>
#include <cassert>

#include <sgm/Pattern.hh>
#include <sgm/Graph_Interface.hh>
#include <sgm/Match_Reporter.hh>


namespace sgm {


//=========================================================================


	 /*! @brief Interface edge match constraint
	  *
	  * A match node constraint describes additional properties that have to be
	  * fulfilled by a given matched edge on a given target.
	  *
	  * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MC_Edge : public Pattern_Interface::Match_Constraint {

	public:

		  //! the source ID of the edge to be constrained
		size_t constrainedFromID;

		  //! the target ID of the edge to be constrained
		size_t constrainedToID;

		  //! construction
		  //! @param constrainedFromID source ID of the edge to be constrained
		  //! @param constrainedToID target ID of the edge to be constrained
		MC_Edge( const size_t constrainedFromID
					, const size_t constrainedToID )
		  :	constrainedFromID(constrainedFromID)
			, constrainedToID(constrainedToID)
		{}

		  //! copy construction
		  //! @param toCopy the MC_Edge object to copy
		MC_Edge( const MC_Edge& toCopy )
		  :	constrainedFromID(toCopy.constrainedFromID)
			, constrainedToID(toCopy.constrainedToID)
		{}

		  //! destruction
		virtual
		~MC_Edge()
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
			return nodeID == constrainedFromID || nodeID == constrainedToID;
		}

		 /*!
		  * Creates a new MC_Edge heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_Edge object that equals this
		  */
		virtual
		MC_Edge*
		clone( void ) const = 0;


		 /*!
		  * Creates a new MC_Edge heap object that equals the current
		  * object but uses the new indices given by old2newIndexMapping.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @param old2newIndexMapping the index mapping to be used for the
		  *        remapping
		  * @param unmatchedIndex an optional specific index that marks
		  *        unmatched nodes within old2newIndexMapping. if this
		  *        constrains one of these nodes, no remapping is done and
		  *        NULL is returned
		  * @return a new allocated MC_Edge object
		  */
		virtual
		MC_Edge*
		remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX ) = 0;


		 /*!
		  * Equality comparison to another match constraint.
		  * @param toCompare the constraint to compare to
		  * @return true if the constraints are equal; false otherwise
		  */
		virtual
		bool
		operator==( const Match_Constraint & toCompare ) const {
			const MC_Edge* pc = dynamic_cast< const MC_Edge* >(&toCompare);
			if (pc != NULL) {
				return this->constrainedFromID == pc->constrainedFromID
						&& this->constrainedToID == pc->constrainedToID;
			} else {
				return false;
			}
		}

	};

//=========================================================================

	  /*! @brief Forbids an edge within a match
	   *
	   * This match constraint defines a given edge is NOT PRESENT within the
	   * target graph.
	   *
	   * Current idea of GML encoding:
	   *
	   * \verbatim
*****************************************************
 constrainNoEdge [
   source   ::= INTEGER
   target   ::= INTEGER
 ]
*****************************************************
	   \endverbatim
	   *
	   * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */
	class MC_NoEdge : public MC_Edge {

	public:

		  //! the index of the source node this constraint is for
		using MC_Edge::constrainedFromID;
		  //! the index of the target node this constraint is for
		using MC_Edge::constrainedToID;

		  /*!
		   * Default construction of the match constraint : undefined values are
		   * initialized with INT_MAX
		   */
		MC_NoEdge()
		 : MC_Edge((size_t)INT_MAX,(size_t)INT_MAX)
		{}

		  /*!
		   * Copy construction of the match constraint
		   */
		MC_NoEdge( const MC_NoEdge& toCopy )
		 : MC_Edge( toCopy )
		{}

		  /*!
		   * Construction of the match constraint
		   *
		   * @param constrainedFromID the index of the source node
		   * @param constrainedToID the index of the target node
		   */
		MC_NoEdge(	const size_t constrainedFromID
					, const size_t constrainedToID )
		 : MC_Edge(constrainedFromID,constrainedToID)
		{}


		  //! destruction
		virtual
		~MC_NoEdge()
		{}

		 /*!
		  * Checks whether or not a match on a given target fullfills the
		  * additional NO-EDGE constraint for the pattern matching.
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
			  // check if no edges existent
			return target.getEdgesBegin(
								match.at(constrainedFromID)
								, match.at(constrainedToID) )
					== target.getEdgesEnd(
								match.at(constrainedFromID)
								, match.at(constrainedToID)
								);
		}

		virtual
		bool
		isConstrainedLabel(	const std::string & label ) const
		{
			  // not label information used
			return false;
		}

		 /*!
		  * Creates a new MC_NoEdge heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_NoEdge object that equals *this
		  */
		virtual
		MC_NoEdge*
		clone( void ) const
		{
			return new MC_NoEdge(*this);
		}


		 /*!
		  * Creates a new MC_NoEdge heap object that equals the current
		  * object but uses the new indices given by old2newIndexMapping.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @param old2newIndexMapping the index mapping to be used for the
		  *        remapping
		  * @param unmatchedIndex an optional specific index that marks
		  *        unmatched nodes within old2newIndexMapping. if this
		  *        constrains one of these nodes, no remapping is done and
		  *        NULL is returned
		  * @return a new allocated MC_NoEdge object
		  */
		virtual
		MC_NoEdge*
		remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
		{
			assert(this->constrainedFromID < old2newIndexMapping.size());
			assert(this->constrainedToID < old2newIndexMapping.size());

			  // check if this edge constrains an unmatched node
			if (old2newIndexMapping.at(this->constrainedFromID) == unmatchedIndex
				|| old2newIndexMapping.at(this->constrainedToID) == unmatchedIndex)
			{
				return NULL;
			}

			  // create copy
			MC_NoEdge* copy = new MC_NoEdge(*this);
			  // do remapping
			copy->constrainedFromID = old2newIndexMapping.at(copy->constrainedFromID);
			copy->constrainedToID = old2newIndexMapping.at(copy->constrainedToID);
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
		operator==( const Match_Constraint& toCompare ) const
		{
			const MC_NoEdge* pc = dynamic_cast< const MC_NoEdge* >(&toCompare);
			if (pc != NULL) {
				bool isEqual = MC_Edge::operator==( toCompare );
				  // return final comparison result
				return isEqual;
			} else {
				return false;
			}
		}



	};

//=========================================================================


	 /*! @brief Constrains edge labels within a match
	  *
	  * A match constraint that restricts the allowed matched labels for a given
	  * matched edge to a set of maintained labels.
	  *
	  * Via the operator it is defined whether the given set of edge labels is
	  * the allowed (op =) or forbidden (op !) set of labels for the edge. It
	  * defaults to allowed (op =).
	  *
	  * NOTE: This constraint is only useful if the according pattern edge shows
	  *       a wildcard label such that it can be matched on any target edge.
	  *
	  * NOTE: if there is NO EDGE between the nodes, the constraint fails!
	  *
	   * Current idea of GML encoding:
	   *
	   * \verbatim
*****************************************************
 constrainEdge [
   source   ::= INTEGER
   target   ::= INTEGER
   op       ::= OPERATOR { = | ! }
   edgeLabels [
     label ::= STRING
   ]
 ]
*****************************************************
	   \endverbatim
	  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MC_EdgeLabel : public MC_Edge {

	public:

		  //! the type that defines a set of labels
		typedef std::set< std::string > LabelSet;

		enum CompareType { ALLOWED, FORBIDDEN };


	public:

		  //! the node ID to be constrained
		using MC_Edge::constrainedFromID;

		  //! the node ID to be constrained
		using MC_Edge::constrainedToID;

		  //! the set of labels this node can be mapped on
		  //! NOTE: if empty, no matching will be allowed !!!
		LabelSet edgeLabels;

		  //! the type of comparison to be applied, i.e. if to ensure that the
		  //! edge label is among the edgeLabels (ALLOWED) or that it is not
		  //! present (FORBIDDEN)
		CompareType compareType;

		  /*!
		   * Default construction of the match constraint : undefined values are
		   * initialized with INT_MAX or empty sets
		   */
		MC_EdgeLabel()
		 : MC_Edge((size_t)INT_MAX,(size_t)INT_MAX)
		   , edgeLabels()
			, compareType(ALLOWED)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedFromID the source node ID to be constrained
		  //! @param constrainedToID the target node ID to be constrained
		MC_EdgeLabel( const size_t constrainedFromID
						, const size_t constrainedToID )
		  :	MC_Edge(constrainedFromID, constrainedToID)
			, edgeLabels()
			, compareType(ALLOWED)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedFromID the source node ID to be constrained
		  //! @param constrainedToID the target node ID to be constrained
		  //! @param edgeLabels the allowed node labels to be matched on
		MC_EdgeLabel( const size_t constrainedFromID
						, const size_t constrainedToID
						, const LabelSet& edgeLabels)
		  :	MC_Edge( constrainedFromID, constrainedToID )
			, edgeLabels( edgeLabels )
			, compareType(ALLOWED)
		{}

		  //! construction where the allowed node labels are empty
		  //! @param constrainedFromID the source node ID to be constrained
		  //! @param constrainedToID the target node ID to be constrained
		  //! @param edgeLabels the allowed node labels to be matched on
		  //! @param compareType the type of comparison to be applied
		MC_EdgeLabel( const size_t constrainedFromID
						, const size_t constrainedToID
						, const LabelSet& edgeLabels
						, const CompareType& compareType )
		  :	MC_Edge( constrainedFromID, constrainedToID )
			, edgeLabels( edgeLabels )
			, compareType( compareType )
		{}

		  //! copy construction
		  //! @param toCopy the MC_EdgeLabel object to copy
		MC_EdgeLabel( const MC_EdgeLabel& toCopy )
		  :	MC_Edge(toCopy)
			, edgeLabels( toCopy.edgeLabels )
			, compareType( toCopy.compareType )
		{}

		  //! destruction
		virtual
		~MC_EdgeLabel()
		{}

		 /*!
		  * Checks whether or not the matched edges between two nodes show
		  * one of the allowed labels or are not showing any of the forbidden
		  * labels.
		  *
		  * NOTE: if there is NO EDGE between the node, the constraint fails!
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
			  // get edge iterator for the edges between source and target
			sgm::Graph_Interface::Edge_iterator curEdge = target.getEdgesBegin(
										match.at(this->constrainedFromID)
										, match.at(this->constrainedToID) );
			sgm::Graph_Interface::Edge_iterator edgeEnd = target.getEdgesEnd(
										match.at(this->constrainedFromID)
										, match.at(this->constrainedToID) );

			  // ensure that there is an edge to constrain ...
			if (curEdge == edgeEnd) {
				return false;
			}

			  // handle comparison according to given comparison type
			switch (compareType) {

			case ALLOWED : { // check if all edge labels are found
				  // empty set of allowed labels -> no label allowed
				if (edgeLabels.empty()) {
					return false;
				}
				  // if wildcard is among the allowed labels -> return true
				if (pattern.getUsedWildcard() != NULL
					&& edgeLabels.find(*(pattern.getUsedWildcard())) != edgeLabels.end() )
				{
					return true;
				}
				  // check if all mapped edges shows an allowed label
				bool allFound = true;
				for( ; allFound && curEdge != edgeEnd; ++curEdge ) {
					allFound = edgeLabels.find( curEdge->getEdgeLabel() )
									!= edgeLabels.end();
				}
				return allFound;
			}
			case FORBIDDEN : { // check that all edge labels are NOT found
				  // empty set of forbidden labels -> all label allowed
				if (edgeLabels.empty()) {
					return true;
				}
				  // if wildcard is among the forbidden pattern -> return false
				if (pattern.getUsedWildcard() != NULL
					&& edgeLabels.find(*(pattern.getUsedWildcard())) != edgeLabels.end() )
				{
					return false;
				}
				  // check if all mapped edges shows no forbidden label
				bool allNotFound = true;
				for( ; allNotFound && curEdge != edgeEnd; ++curEdge ) {
					allNotFound = edgeLabels.find( curEdge->getEdgeLabel() )
									== edgeLabels.end();
				}
				return allNotFound;
			}
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
			  // check if among the edge labels
			for (LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end(); ++l)
				if (l->compare(label)==0)
					return true;
			  // not found
			return false;
		}

		 /*!
		  * Creates a new MC_EdgeLabel heap object that equals the current
		  * object.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @return a new allocated MC_EdgeLabel object that equals *this
		  */
		virtual
		MC_EdgeLabel*
		clone( void ) const
		{
			return new MC_EdgeLabel(*this);
		}


		 /*!
		  * Creates a new MC_EdgeLabel heap object that equals the current
		  * object but uses the new indices given by old2newIndexMapping.
		  * NOTE: YOU have to delete it later on! There is no garbage
		  *       collection!
		  * @param old2newIndexMapping the index mapping to be used for the
		  *        remapping
		  * @param unmatchedIndex an optional specific index that marks
		  *        unmatched nodes within old2newIndexMapping. if this
		  *        constrains one of these nodes, no remapping is done and
		  *        NULL is returned
		  * @return a new allocated MC_EdgeLabel object
		  */
		virtual
		MC_EdgeLabel*
		remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
		{
			assert(this->constrainedFromID < old2newIndexMapping.size());
			assert(this->constrainedToID < old2newIndexMapping.size());

			  // check if this edge constrains an unmatched node
			if (old2newIndexMapping.at(this->constrainedFromID) == unmatchedIndex
				|| old2newIndexMapping.at(this->constrainedToID) == unmatchedIndex)
			{
				return NULL;
			}

			  // create copy
			MC_EdgeLabel* copy = new MC_EdgeLabel(*this);
			  // do remapping
			copy->constrainedFromID = old2newIndexMapping.at(copy->constrainedFromID);
			copy->constrainedToID = old2newIndexMapping.at(copy->constrainedToID);
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
		operator==( const Match_Constraint& toCompare ) const
		{
			const MC_EdgeLabel* pc = dynamic_cast< const MC_EdgeLabel* >(&toCompare);
			if (pc != NULL) {
				bool isEqual =
					MC_Edge::operator==( toCompare )
					&& this->edgeLabels.size() == pc->edgeLabels.size()
					&& this->compareType == pc->compareType
					;
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

} // namespace sgm


#endif /* SGM_MC_EDGE_HH_ */

