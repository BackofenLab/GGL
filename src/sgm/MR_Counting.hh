#ifndef SGM_MR_COUNTING_HH_
#define SGM_MR_COUNTING_HH_

#include "sgm/Match_Reporter.hh"


namespace sgm {
	
	  /*! @brief Counts matches
	   *
	   *  A sgm::Match_Reporter implementation that only counts the number of
	   *  reported matches.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class MR_Counting : public Match_Reporter {
		
	protected:
		
		  //! the number of matches reported so far
		size_t numberOfMatches;
		
	public:
		
		  //! Construction
		MR_Counting( );
		
		virtual 
		~MR_Counting();
		
		  //! Increases the internal counter, nothing else.
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
		
		  //! Resets the number of matches reported so far to 0.
		void
		resetHits ( void );
		
		  //! Access to the number of reported matches so far
		  //! @return number of matches
		size_t
		getHits( void ) const;
	
	};
	
} // namespace sgm

 // IMPLEMENTATION OF MR_Counting

namespace sgm {

	inline
	MR_Counting
	::MR_Counting( )
	 : numberOfMatches(0)
	{}
	
	inline
	MR_Counting
	::~MR_Counting( )
	{}

	inline
	void
	MR_Counting
	::reportHit (	const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match )
	{
		numberOfMatches++;
	}

	inline
	void
	MR_Counting
	::resetHits ( void )
	{
		numberOfMatches = 0;
	}

	inline
	size_t
	MR_Counting
	::getHits ( void ) const
	{
		return numberOfMatches;
	}

}

#endif /*MR_COUNTING_HH_*/
