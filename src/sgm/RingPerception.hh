
#ifndef SGM_RINGPERCEPTION_HH_
#define SGM_RINGPERCEPTION_HH_

#include "sgm/Graph_Interface.hh"
#include "sgm/RingReporter.hh"

namespace sgm {

	  /*! @brief Interface ring enumeration
	   *
	   * Generic interface for ring perception algorithms that allow for the
	   * enumeration of all rings within a graph.
	   *
	   */
	class RingPerception {
	public:

		  //! construction
		RingPerception() {}
		  //! destruction
		virtual ~RingPerception() {}

		  /*!
		   * Finds rings within the given graph and reports each found ring to
		   * the assigned reporter.
		   *
		   * @param graph the graph to be analyzed
		   * @param reporter the RingReporter to report all found rings to
		   * @param maxRingSize the maximal size of rings to report
		   * @return the number of rings found and reported
		   */
		virtual
		size_t
		findRings( const Graph_Interface & graph
					, RingReporter & reporter
					, const size_t maxRingSize = std::numeric_limits<size_t>::max() ) = 0;

	};

} // namespace sgm

#endif /* SGM_RINGPERCEPTION_HH_ */
