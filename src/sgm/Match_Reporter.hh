#ifndef SGM_MATCH_REPORTER_HH_
#define SGM_MATCH_REPORTER_HH_

#include "sgm/Match.hh"
#include "sgm/Graph_Interface.hh"
#include "sgm/Pattern.hh"

namespace sgm {

	  /*! @brief Interface graph match reporting
	   *
	   *  An interface description of the class used by sgm::SubGraphMatching objects
	   *  to report a match of the pattern in the target graph.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Match_Reporter {
		
	public:
		
		virtual
		~Match_Reporter()
		{}
		
		  //! Reports a match. The match is encoded using a vector. The length 
		  //! of the vector corresponds to the number of vertices in the pattern
		  //! and position i encodes the matched position of pattern node i in
		  //! the target graph.
		  //! @param pattern the pattern graph that was searched for
		  //! @param target the graph the pattern was found within
		  //! @param match contains the indices of the matched pattern nodes in
		  //! the target graph. match[i] corresponds to the mapping of the ith
		  //! vertex in the pattern graph.
		virtual
		void
		reportHit (	const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match ) = 0;
		
	};

} // namespace sgm


#endif /*MATCH_REPORTER_HH_*/
