#ifndef SGM_MR_STREAM_HH_
#define SGM_MR_STREAM_HH_

#include "sgm/Match_Reporter.hh"

#include <iostream>
#include <iomanip>

namespace sgm {
	
	  /*! @brief Writes each match to an output stream
	   *
	   *  A sgm::Match_Reporter implementation that writes the matched pairs of
	   *  pattern and target graph nodes to stream.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class MR_stream : public Match_Reporter {
		
	protected:
		
		  //! the stream to write to
		std::ostream & out;
		
		  //! the number of matches reported so far
		size_t numberOfMatches;
		
	public:
		
		  //! Construction
		MR_stream( std::ostream& out );
		
		virtual 
		~MR_stream();
		
		  //! Writes a match to stream.
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
	
	};
	
} // namespace sgm

 // IMPLEMENTATION OF MR_stream

namespace sgm {

	inline
	MR_stream
	::MR_stream( std::ostream& out_ )
	 : out(out_), numberOfMatches(0)
	{}
	
	inline
	MR_stream
	::~MR_stream( )
	{}

	inline
	void
	MR_stream
	::reportHit (	const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match )
	{
		numberOfMatches++;
		
		out <<std::setw(6) <<numberOfMatches <<" :";
		for (size_t i=0; i<match.size(); ++i) {
			out <<" (" <<i <<"," <<match[i] <<")";
		}
		out <<std::endl;
	}

	inline
	void
	MR_stream
	::resetHits ( void )
	{
		numberOfMatches = 0;
	}

}

#endif /*MR_STREAM_HH_*/
