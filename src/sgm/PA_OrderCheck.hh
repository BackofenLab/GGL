#ifndef SGM_PA_ORDERCHECK_HH_
#define SGM_PA_ORDERCHECK_HH_


#include <algorithm>

#include "sgm/Pattern_Automorphism.hh"
#include "sgm/MR_Storing.hh"
#include "sgm/GM_vf2.hh"


namespace sgm {


	  /*! @brief Order based automorphism exclusion
	   *
	   * Performs an order check (e.g. "<") on match positions to identify 
	   * symmetric matches.
	   */
	class PA_OrderCheck : public Pattern_Automorphism {
	public:
		
		  //! pair of indices for order checking
		typedef std::pair<	Graph_Interface::IndexType
							, Graph_Interface::IndexType > IndexPair;
		  //! list of Match positions to do an order check on
		typedef std::vector< IndexPair > CheckList;
		
	protected:
		
		  //! the source graph this automorphism description is defined for
		const Pattern_Interface* pattern;
		
		  //! the list of Match positions to do an order check on
		const CheckList checkList;

	public:
		
		 /*! Construction of a graph symmetry handler.
		  * @param pattern the graph this Pattern_Automorphism handles
		  * @param checkList the list of order checks to apply to break all
		  *        symmetries in the graph
		  */
		PA_OrderCheck(	const Pattern_Interface& pattern
						, const CheckList& checkList
						 );
		
		 /*! Copy construction of a graph symmetry handler.
		  * @param toCopy the object to make this a copy of
		  */
		PA_OrderCheck(	const PA_OrderCheck & toCopy );

		 //! Destruction
		virtual
		~PA_OrderCheck();
		
		 /*!
		  * Checks whether or not the given match is symmetric to another match
		  * that is possible based on a list of order checks on the match.
		  * 
		  * @param graph the graph the match belongs to
		  * @param match the match of the graph that has to be checked
		  * @return true if the match is a symmetric one, false otherwise
		  */
		bool
		isSymmetryMatch(	const Pattern_Interface& graph
							, const Match& match ) const;
		
		 /*! Access to the internally used list of order checks applied in the
		  * symmetry check.
		  * @return the used CheckList
		  */
		const CheckList&
		getCheckList(void) const;
		
		 /*! Access to the graph this Pattern_Automorphism belongs to.
		  * @return the source graph
		  */
		const Pattern_Interface&
		getPattern(void) const;
		
		 /*! Creates a PA_OrderCheck object that handles all symmetries of a
		  * given graph. The template class specifies the GraphMatcher to
		  * be used.
		  *
		  * NOTE: PA_OrderCheck instances are only valid for patterns WITHOUT
		  * additional matching constraints. Thus, the method will return no
		  * checks for patterns including constraints!
		  *
		  * @param pattern the pattern graph of interest
		  * @return the Pattern_Automorphism breaking all symmetries of the graph
		  */
		template < class GRAPHMATCHER >
		static
		PA_OrderCheck
		getGraphAutomorphismT( const Pattern_Interface& pattern );

		 /*! Creates a PA_OrderCheck object that handles all symmetries of a
		  * given graph. Calls getGraphAutomorphismT< GM_vf2 >.
		  *
		  * NOTE: PA_OrderCheck instances are only valid for patterns WITHOUT
		  * additional matching constraints. Thus, the method will return no
		  * checks for patterns including constraints!
		  *
		  * @param pattern the pattern graph of interest
		  * @return the Pattern_Automorphism breaking all symmetries of the graph
		  */
		static
		PA_OrderCheck
		getGraphAutomorphism( const Pattern_Interface& pattern );
	

		 /*!
		  * Creates a new allocated copy of this object on heap. NOTE, the
		  * returned copy has to be deleted by the calling method.
		  *
		  * @return a copy of this
		  */
		virtual
		PA_OrderCheck *
		clone( void ) const;

	};
	

} // namespace sgm


#include <ostream>

/*!
 * Prints the list of order checks of a PA_OrderCheck object to stream.
 * 
 * @param out the stream to write to
 * @param pa the PA_OrderCheck object to write
 * @return the modified stream
 */
std::ostream&
operator <<( std::ostream& out, const sgm::PA_OrderCheck& pa );


// include implementation
#include "sgm/PA_OrderCheck.icc"

#endif /*SGM_PA_ORDERCHECK_HH_*/
