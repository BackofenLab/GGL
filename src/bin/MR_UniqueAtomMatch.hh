
#ifndef MR_UNIQUEATOMMATCH_H_
#define MR_UNIQUEATOMMATCH_H_


#include <sgm/Match_Reporter.hh>

/**
 * If the target graph is an instance of
 *   ggl::chem::Molecule_Graph_V
 * or such an instance wrapped by
 *   ggl::chem::Molecule_Graph_noClass
 * than only matches are forwarded that match on different atoms of the
 * underlying molecule(s).
 *
 */
class MR_UniqueAtomMatch: public sgm::Match_Reporter
{
public:

	/**
	 * Construction
	 *
	 * @param forwardReporter the MR to report validated matches to
	 */
	MR_UniqueAtomMatch( sgm::Match_Reporter & forwardReporter );

	/**
	 * destruction
	 */
	virtual ~MR_UniqueAtomMatch();


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

protected:

	//! the MR to report validated matches to
	sgm::Match_Reporter & forwardReporter;

};

#endif /* MRUNIQUEATOMMATCH_H_ */
