
#ifndef SGM_RINGREPORTER_HH_
#define SGM_RINGREPORTER_HH_

#include <set>
#include <list>

#include "sgm/Graph_Interface.hh"

namespace sgm {

	  /*! @brief Interface to report found rings
	   *
	   * Abstract class that describes the interface to report rings in graphs.
	   * Thus the interface is used by sgm::RingPerception to report found
	   * rings.
	   *
	   * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class RingReporter {

	public:

		  //! type of a set of node indices that form a ring
		typedef std::set<size_t> RingNodes;
		  //! type of a list of node indices that form a ring
		typedef std::list<size_t> RingList;

	public:
		  //! construction
		RingReporter();
		  //! destruction
		virtual ~RingReporter();

		  /*!
		   * Is called to report a ring.
		   * @param graph the graph that contains the ring
		   * @param ringList the ring to report
		   */
		virtual
		void
		reportRing( const Graph_Interface& graph, const RingList & ringList ) = 0;

	public:

		 /*!
		  * Utility class that converts a RingList, i.e. a list of nodes when
		  * traversing a ring within a graph, into the set of nodes that form
		  * the ring.
		  *
		  * @param graph the graph that contains the ring
		  * @param ringList the ring list to convert
		  * @return the list of nodes traversing the ring
		  */
		static
		RingNodes
		toRing( const Graph_Interface& graph, const RingList& ringList );
	};

} // namespace sgm

#endif /* SGM_RINGREPORTER_HH_ */
