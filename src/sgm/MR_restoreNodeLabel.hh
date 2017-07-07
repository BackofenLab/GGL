#ifndef MR_RESTORENODELABEL_HH_
#define MR_RESTORENODELABEL_HH_

namespace sgm
{



//////////////////////////////////////////////////////////////////////////

#include <sgm/Match_Reporter.hh>

/**
 * Wrapper class to restore the node label information of the target graph that
 * was altered using sgm::Graph_NodeLabelPrefix, which is then forwarded to
 * a given Match_Reporter instance.
 *
 */
class MR_restoreNodeLabel : public sgm::Match_Reporter {

protected:

	//! the Match_Reporter to forward the match information to with full node
	//! label information for the target graph
	sgm::Match_Reporter & forwardMR;

public:

	//! construction
	//! @param forwardMR the Match_Reporter to forward the match information
	//!            to with full node label information for the target graph
	MR_restoreNodeLabel( sgm::Match_Reporter & forwardMR );

	//! destruction
	virtual
	~MR_restoreNodeLabel();

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
	reportHit (	const sgm::Pattern_Interface & pattern,
				const sgm::Graph_Interface & target,
				const sgm::Match & match );


};



} /* namespace sgm */


#include "sgm/MR_restoreNodeLabel.icc"


#endif /* MR_RESTORENODELABEL_HH_ */
