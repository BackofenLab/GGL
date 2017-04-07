#ifndef SGM_PATTERN_AUTOMORPHISM_HH_
#define SGM_PATTERN_AUTOMORPHISM_HH_

#include "sgm/Match.hh"
#include "sgm/Pattern.hh"


namespace sgm {

	/*! @brief Automorphism description for patterns
	 *
	 * An encoding of what matches of a pattern are symmetric to other.
	 * 
	 * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */
	class Pattern_Automorphism {
		
	public:
		
		virtual
		~Pattern_Automorphism()
		{}
		
		 /*!
		  * Checks whether or not the given match is symmetric to another match
		  * that is possible.
		  *
		  * @param graph the graph the match belongs to
		  * @param match the match of the graph that has to be checked 
		  * @return true if the match is a symmetric one, false otherwise
		  */
		virtual
		bool
		isSymmetryMatch(	const Pattern_Interface& graph
							, const Match& match ) const = 0;
		
		 /*!
		  * Creates a new allocated copy of this object on heap. NOTE, the
		  * returned copy has to be deleted by the calling method.
		  *
		  * @return a copy of this
		  */
		virtual
		Pattern_Automorphism *
		clone( void ) const = 0;

	};

} // namespace sgm

#endif /* SGM_PATTERN_AUTOMORPHISM_HH_*/
