
#ifndef SGM_MATCH_HH_
#define SGM_MATCH_HH_

#include "sgm/Graph_Interface.hh"

namespace sgm {

	//! encoding of a match, i.e. for each node index of the pattern the
	//! matched node within the target graph is stored
	typedef std::vector<Graph_Interface::IndexType> Match;



}  // namespace sgm

#endif /* SGM_MATCH_HH_ */
