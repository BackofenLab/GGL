#ifndef SGM_PATTERN_HH_
#define SGM_PATTERN_HH_

#include "sgm/Graph_Interface.hh"
#include "sgm/Match.hh"
#include <climits>


namespace sgm {


	  /*! @brief Pattern description to be matched
	   *
	   *  Abstract pattern description to be used in graph matching algorithms.
	   *
	   *  @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Pattern_Interface {

	public:

	//=========================================================================


		 /*!
		  * A match constraint describes additional properties that have to be
		  * fullfilled by a given match on a given target.
		  *
		  */
		class Match_Constraint {

		public:

			  //! construction
			Match_Constraint(  )
			{}

			  //! destruction
			virtual
			~Match_Constraint()
			{}

			 /*!
			  * Checks whether or not a match on a given target fullfills the
			  * additional matching constraint.
			  *
			  * @param pattern the pattern graph that was matched
			  * @param target the target graph the pattern was matched on
			  * @param match the match information for the pattern
			  *        on the target graph
			  * @return true if the match is valid; false if the constraint is
			  *         violated
			  */
			virtual
			bool
			isValidMatch(	const Pattern_Interface & pattern,
							const Graph_Interface & target,
							const Match & match ) const = 0;

			 /*!
			  * Checks whether or not a given label is part of the constraint
			  * information. This check is needed by some parsers to verify the
			  * wildcard definition.
			  *
			  * @param label the label of interest
			  * @return true if the label is part of the constraint; false otherwise
			  */
			virtual
			bool
			isConstrainedLabel(	const std::string & label ) const = 0;

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
			isConstraining(	const size_t nodeID ) const = 0;


			 /*!
			  * Creates a new Match_Constraint heap object that equals the current
			  * object.
			  * NOTE: YOU have to delete it later on! There is no garbage
			  *       collection!
			  * @return a new allocated Match_Constraint object that equals this
			  */
			virtual
			Match_Constraint*
			clone( void ) const = 0;

			 /*!
			  * Creates a new Match_Constraint heap object that equals the current
			  * object but uses the new indices given by old2newIndexMapping.
			  * NOTE: YOU have to delete it later on! There is no garbage
			  *       collection!
			  * @param old2newIndexMapping the index mapping to be used for the
			  *        remapping
			  * @param unmatchedIndex an optional specific index that marks
			  *        unmatched nodes within old2newIndexMapping. if this
			  *        constrains one of these nodes, no remapping is done and
			  *        NULL is returned
			  * @return a new allocated Match_Constraint object or NULL if the
			  *        constraint was covering unmatched nodes
			  */
			virtual
			Match_Constraint*
			remap( const Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX ) = 0;

			 /*!
			  * Equality comparison to another match constraint.
			  * @param toCompare the constraint to compare to
			  * @return true if the constraints are equal; false otherwise
			  */
			virtual
			bool
			operator==( const Match_Constraint& toCompare ) const = 0;

		};

	//=========================================================================


		  //! container of matching constraints
		typedef std::vector< Match_Constraint* > ConstraintVec;


	public:


		  //! destruction
		virtual ~Pattern_Interface();

		  //! Equality comparison
		  //! @param toCompare the Pattern to compare to
		  //! @return true if both describe the same pattern
		virtual
		bool
		operator==(const Pattern_Interface& toCompare ) const;

		  //! Inequality comparison
		  //! @param toCompare the Pattern to compare to
		  //! @return true if both describe the different patterns
		virtual
		bool
		operator!=(const Pattern_Interface& toCompare ) const;


		  //! Access to the pattern graph to be matched
		  //! @return the pattern graph
		virtual
		const Graph_Interface &
		getPatternGraph( void ) const = 0;

		  //! Access to the matching constraints that have to be fulfilled by
		  //! a match of this pattern graph
		  //! @return the matching constraints to validate
		virtual
		const ConstraintVec&
		getConstraints( void ) const = 0;

		  //! Access to the wildcard to be used when matching this pattern onto
		  //! some other graph.
		  //! @return the wildcard string to be used for edge and node labels,
		  //!         or NULL if no wildcard should be applied
		virtual
		const std::string *
		getUsedWildcard( void ) const = 0;

	};

} // namespace sgm



namespace sgm {


	  //! A wrapper class around a given Graph_interface that can serve as a
	  //! pattern.
	  //!
	  //! @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	class Pattern : public Pattern_Interface {

	protected:

		  //! the graph to be represented as a pattern
		const Graph_Interface * graph;

		  //! the additional match constraints to be fulfilled by each match
		ConstraintVec matchConstraints;

		  //! the wildcard string to be used for matching
		const std::string * usedWildcard;

	public:

		  //! constructs a pattern without additional constraints
		  //! @param graph the graph to represent
		Pattern( const Graph_Interface& graph );

		  //! constructs a pattern without additional constraints
		  //! @param graph the graph to represent
		  //! @param wildcardToUse the wildcard to use for matching
		Pattern( const Graph_Interface& graph
				, const std::string &  wildcardToUse );

		  //! constructs a pattern with additional constraints
		  //! @param graph the graph to represent
		  //! @param matchConstraints the additional constraints to be fulfilled
		Pattern( const Graph_Interface& graph
				, const ConstraintVec & matchConstraints );

		  //! constructs a pattern with additional constraints
		  //! @param graph the graph to represent
		  //! @param matchConstraints the additional constraints to be fulfilled
		  //! @param wildcardToUse the wildcard to use for matching
		Pattern( const Graph_Interface& graph
				, const ConstraintVec & matchConstraints
				, const std::string & wildcardToUse );

		  //! copy construction
		  //! @param toCopy the pattern to make this object a copy of
		Pattern( const Pattern& toCopy );

		  //! destruction
		virtual ~Pattern();


		  //! Access to the pattern graph to be matched
		  //! @return the pattern graph
		virtual
		const Graph_Interface &
		getPatternGraph( void ) const;

		  //! Access to the matching constraints that have to be fulfilled by
		  //! a match of this pattern graph
		  //! @return the matching constraints to validate
		virtual
		const ConstraintVec&
		getConstraints( void ) const;

		  //! Access to the wildcard to be used when matching this pattern onto
		  //! some other graph.
		  //! @return the wildcard string to be used for edge and node labels,
		  //!         or NULL if no wildcard should be applied
		virtual
		const std::string *
		getUsedWildcard( void ) const;

		  //! assignment operator
		  //! @param toCopy the instance to make this a copy of
		  //! @return the altered object (*this)
		sgm::Pattern&
		operator=(const sgm::Pattern& toCopy );

	public:

		  //! Does a deep copy of the given constraint vector
		  //! @param toCopy the ConstraintVec instance to copy
		  //! @return the deep copy of the given vector of MatchConstraints
		static
		ConstraintVec
		copyConstraintVec( const ConstraintVec & toCopy );

	};

} // namespace sgm



#include <iostream>
#include <iomanip>


  /*! Prints a Pattern instance to stream. For each node its label and
   * the adjancent nodes including the edge label is printed.
   * 
   * @param out the stream to write to
   * @param g the graph to write
   * @return the modified out stream
   */
std::ostream&
operator <<( std::ostream & out, const sgm::Pattern& g );


#endif /*SGM_PATTERN_HH_*/
