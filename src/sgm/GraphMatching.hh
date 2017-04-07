#ifndef SGM_GRAPHMATCHING_HH_
#define SGM_GRAPHMATCHING_HH_

#include "sgm/Graph_Interface.hh"
#include "sgm/Pattern.hh"
#include "sgm/Match_Reporter.hh"

namespace sgm {


	  /*! @brief Interface graph isomorphism
	   *
	   *  This class defines the interface of graph matching algorithms of
	   *  the sgm library, i.e. algorithms solving the graph-isomorphism
	   *  problem.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
class GraphMatching {
		
	public:
		
		virtual
		~GraphMatching()
		{}

		  //! Performs exact graph matching to find maxHits occurrences of the
		  //! pattern graph within the target graph. Each hit is reported to
		  //! the Match_Reporter object.
		  //! @param pattern the pattern graph to search for
		  //! @param target the graph to search the pattern within
		  //! @param reporter all hits are reported to that object
		  //! @param maxHits the maximal number of hits to find
		  //! @return the number of exact matches found
		virtual
		size_t
		findMatches (	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						Match_Reporter & reporter,
						const size_t maxHits ) = 0;

		  //! Performs exact graph matching to find maxHits occurrences of the
		  //! pattern graphs within the target graph. Each hit is reported to
		  //! the Match_Reporter object.
		  //! NOTE, the first maxHits matches of the pattern graphs according
		  //! to their order in the patterns container are reported, i.e. first
		  //! all occurrences of the first pattern are identified. If this does
		  //! not exceed the maxHits limit, the next pattern is matched and so
		  //! on until either no pattern is left or the maxHits limit is
		  //! exceeded.
		  //! @param patterns the container of the pattern graphs to search for
		  //! @param target the graph to search the pattern within
		  //! @param reporter all hits are reported to that object
		  //! @param maxHits the maximal number of hits to find
		  //! @return the number of exact matches found
		virtual
		size_t
		findMatches (	const std::vector< const Pattern_Interface*> & patterns,
						const Graph_Interface & target,
						Match_Reporter & reporter,
						const size_t maxHits ) = 0;

		  //! Performs exact graph matching to find maxHits occurrences of the
		  //! pattern graphs within the target graph. Each hit is reported to
		  //! the Match_Reporter object.
		  //! NOTE, the first maxHits matches of the pattern graphs according
		  //! to their order in the patterns container are reported, i.e. first
		  //! all occurrences of the first pattern are identified. If this does
		  //! not exceed the maxHits limit, the next pattern is matched and so
		  //! on until either no pattern is left or the maxHits limit is
		  //! exceeded.
		  //! @param patterns the container of the pattern graphs to search for
		  //! @param target the graph to search the pattern within
		  //! @param reporters each hit is reported to the corresponding object,
		  //!        the container has to have the same length as patterns
		  //! @param maxHits the maximal number of hits to find
		  //! @return the number of exact matches found
		virtual
		size_t
		findMatches (	const std::vector< const Pattern_Interface*> & patterns,
						const Graph_Interface & target,
						std::vector<Match_Reporter*> & reporters,
						const size_t maxHits ) = 0;
		
		
	};

} // namespace sgm

#endif /*GRAPHMATCHING_HH_*/

