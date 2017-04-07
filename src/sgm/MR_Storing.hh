#ifndef SGM_MR_STORING_HH_
#define SGM_MR_STORING_HH_

#include "sgm/Match_Reporter.hh"


namespace sgm {
	
	  /*! @brief Stores each match in an STL container
	   *
	   *  A sgm::Match_Reporter implementation that stores each match mapping
	   *  within a provided STL container.
	   *
	   *  STL_PUSHBACK_CONTAINER should allow for "push_back(Match)".
	   *
	   *  @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	template <class STL_PUSHBACK_CONTAINER>
	class MR_StoringT: public Match_Reporter {
		
	public:
		
		typedef STL_PUSHBACK_CONTAINER Storage;
		
	protected:
		
		  //! where to store the matches in using its "push_back" method
		Storage & storage;
		
	public:
		
		  //! Construction
		  //! @param storage the STL_PUSHBACK_CONTAINER to write to
		MR_StoringT( STL_PUSHBACK_CONTAINER & storage );
		
		virtual 
		~MR_StoringT();
		
		  //! Adds the match to the STL container without further processing.
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
	
	//! Template wrapper for default values
	typedef MR_StoringT< std::vector< sgm::Match > > MR_Storing;
	
} // namespace sgm

 // IMPLEMENTATION OF MR_StoringT

namespace sgm {

	template <class STL_PUSHBACK_CONTAINER>
	inline
	MR_StoringT<STL_PUSHBACK_CONTAINER>
	::MR_StoringT( STL_PUSHBACK_CONTAINER & storage_ )
	 :	storage(storage_)
	{}
	
	template <class STL_PUSHBACK_CONTAINER>
	inline
	MR_StoringT<STL_PUSHBACK_CONTAINER>
	::~MR_StoringT( )
	{}

	template <class STL_PUSHBACK_CONTAINER>
	inline
	void
	MR_StoringT<STL_PUSHBACK_CONTAINER>
	::reportHit (	const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match )
	{
		  // store match
		storage.push_back(match);
	}

}


#include <set>

namespace sgm {
	
	  //! A sgm::Match_Reporter implementation that stores each match mapping
	  //! within a provided STL container.
	  //! 
	  //! STL_INSERT_CONTAINER should allow for "push_back(Match)".
	  //!
	  //! @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	template <class STL_INSERT_CONTAINER>
	class MR_StoringInsertT: public Match_Reporter {
		
	public:
		
		typedef STL_INSERT_CONTAINER Storage;
		
	protected:
		
		  //! where to store the matches in using its "push_back" method
		Storage & storage;
		
	public:
		
		  //! Construction
		  //! @param storage the STL_INSERT_CONTAINER to write to
		MR_StoringInsertT( STL_INSERT_CONTAINER & storage );
		
		virtual 
		~MR_StoringInsertT();
		
		  //! Adds the match to the STL container without further processing.
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
	
	//! Template wrapper for default values
	typedef MR_StoringInsertT< std::set< sgm::Match > > MR_StoringInsert;
	
} // namespace sgm

 // IMPLEMENTATION OF MR_StoringInsertT

namespace sgm {

	template <class STL_INSERT_CONTAINER>
	inline
	MR_StoringInsertT<STL_INSERT_CONTAINER>
	::MR_StoringInsertT( STL_INSERT_CONTAINER & storage_ )
	 :	storage(storage_)
	{}
	
	template <class STL_INSERT_CONTAINER>
	inline
	MR_StoringInsertT<STL_INSERT_CONTAINER>
	::~MR_StoringInsertT( )
	{}

	template <class STL_INSERT_CONTAINER>
	inline
	void
	MR_StoringInsertT<STL_INSERT_CONTAINER>
	::reportHit (	const Pattern_Interface & pattern,
					const Graph_Interface & target,
					const Match & match )
	{
		  // store match
		storage.insert(match);
	}

}

#endif /* SGM_MR_STORING_HH_ */

