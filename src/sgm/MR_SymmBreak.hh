#ifndef SGM_MR_SYMMBREAK_HH_
#define SGM_MR_SYMMBREAK_HH_

#include "sgm/Match_Reporter.hh"
#include "sgm/Pattern_Automorphism.hh"


namespace sgm {
	
	  /*! @brief Symmetry breaking among matches
	   *
	   * An sgm::Match_Reporter implementation that wraps another
	   * Match_Reporter and utilizes an sgm::Pattern_Automorphism to determine
	   * if a reported match is a symmetric solution or not. If not the match 
	   * is forwarded to the provided second reporter.
	   * 
	   * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class MR_SymmBreak : public Match_Reporter {
		
	protected:
		
		  //! the Pattern_Automorphism used to check for symmetries
		const Pattern_Automorphism * symmCheck;
		  //! the reporter to forward all non-symmetric matches to
		Match_Reporter* forward;
		
	public:
		
		  /*! Construction
		   * 
		   * @param symmCheck the automorphism checker to apply to the pattern
		   *                  and match reported.
		   * @param forward the reporter to forward all non-symmetric matches to
		   */
		MR_SymmBreak(	const Pattern_Automorphism& symmCheck
						, Match_Reporter& forward );
		
		  /*! Copy construction
		   *
		   * @param toCopy the object to make this a copy of
		   */
		MR_SymmBreak(	const MR_SymmBreak & toCopy );

		  //! Destruction
		virtual 
		~MR_SymmBreak();
		
		MR_SymmBreak &
		operator=( const MR_SymmBreak& toCopy );

		  //! Checks if the match is symmetric and if not forwards the match
		  //! to the other provided reporter. 
		  //! @param pattern the pattern graph that was searched for
		  //! @param target the graph the pattern was found within
		  //! @param match contains the indices of the matched pattern nodes in
		  //! the target graph. match[i] corresponds to the mapping of the ith
		  //! vertex in the pattern graph.
		virtual
		void
		reportHit ( const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match );
		
	
	};
	
} // namespace sgm

 // IMPLEMENTATION OF MR_SymmBreak

namespace sgm {

	inline
	MR_SymmBreak
	::MR_SymmBreak(	const Pattern_Automorphism& symmCheck_
					, Match_Reporter& forward_ )
	 :	symmCheck(symmCheck_.clone())
	 	, forward(&forward_)
	{}

	inline
	MR_SymmBreak
	::MR_SymmBreak(	const MR_SymmBreak& toCopy )
	 :	symmCheck(toCopy.symmCheck->clone())
	 	, forward(toCopy.forward)
	{}
	
	inline
	MR_SymmBreak
	::~MR_SymmBreak( )
	{
		if (symmCheck != NULL) delete symmCheck;
		symmCheck = NULL;
	}

	inline
	MR_SymmBreak &
	MR_SymmBreak
	::operator=( const MR_SymmBreak& toCopy )
	{
		 // deep copy
		if (symmCheck != NULL)
			delete symmCheck;
		symmCheck = toCopy.symmCheck->clone();
		 // flat copy
	 	forward = toCopy.forward;
	 	 // access to the changed *this object
	 	return *this;
	}

	inline
	void
	MR_SymmBreak
	::reportHit (	const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match )
	{
		  // check if symmetric
		if (!symmCheck->isSymmetryMatch( pattern, match )) {
			  // non-symmetric hits are forwarded
			forward->reportHit( pattern, target, match );
		}
	}

}

#endif /*MR_SYMMBREAK_HH_*/
