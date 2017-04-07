
#ifndef SGM_RP_HANSER96_HH_
#define SGM_RP_HANSER96_HH_


#include <boost/graph/adjacency_list.hpp>

#include "sgm/RingPerception.hh"

#include <limits>
#include <set>


namespace sgm {

	 /*! @brief Ring enumeration ala Hanser et al.
	  *
	  * Implementation of the exhaustive ring perception algorithm by
	  * Hanser et al. (1996)
	  *
	  *   A New Algorithm for Exhaustive Ring Perception in a Molecular Graph
	  *   T. Hanser, P. Jauffret, and G. Kaufmann
	  *   J. Chem. Inf. Comput. Sci., 1996, 36, 1146-1152
	  *
	  * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class RP_Hanser96 : public RingPerception {
	public:

		typedef std::set< std::pair<size_t,size_t> > BondSet;

		  //! construction
		RP_Hanser96();

		  //! destruction
		virtual ~RP_Hanser96();


		  /*!
		   * Finds ALL rings within the given graph and reports each found
		   * ring to the assigned reporter. The enumeration can be restricted
		   * to rings up to a given ring size, ie. number of nodes per ring.
		   *
		   * @param graph the graph to be analyzed
		   * @param reporter the RingReporter to report all found rings to
		   * @param maxRingSize the maximal size of rings to report
		   * @return the number of all rings within the gaph
		   */
		size_t
		findRings( const Graph_Interface & graph,
					RingReporter & reporter,
					const size_t maxRingSize = std::numeric_limits<size_t>::max());

		  /*!
		   * Identifies all bonds participating in rings and stores the
		   * according vertex index pairs in the provided container.
		   *
		   * @param graph the graph to search for ring bonds
		   * @param ringBonds OUT : the container to fill with the ring bonds.
		   *   NOTE: the container is cleared at the beginning of the call.
		   * @param maxRingSize the maximal size of rings to consider
		   *
		   * @return the number of ring bonds within the graph
		   */
		size_t
		findRingBonds( const Graph_Interface & graph
					, BondSet & ringBonds
					, const size_t maxRingSize = std::numeric_limits<size_t>::max() );


	protected:

		  //! information on the current ring part
		struct RingInfo {
		public:
			RingReporter::RingNodes nodes;
			RingReporter::RingList path;
		};


		  //! The properties available for the nodes of the P-graph
		typedef	boost::property<	boost::vertex_index_t, size_t >
						Graph_NodeProperties;

		  //! The properties available for the edges of the P-graph
		typedef	boost::property<	boost::edge_name_t, RingInfo >
						Graph_EdgeProperties;

		  //! the P-graph type that is used to detect all rings
		typedef boost::adjacency_list<
							boost::multisetS,      			// store edges
							boost::vecS,       				// store vertices
							boost::undirectedS,				// is an undirected graph
							Graph_NodeProperties,  			// index information
							Graph_EdgeProperties   			// ring information
						>
							P_Graph;

	protected:

		  //! the P-graph to be compressed during the ring perception
		P_Graph pGraph;

		  //! the index access for pGraph
		boost::property_map<P_Graph, boost::vertex_index_t>::type pGraphIndex;
		  //! the edge label access for pGraph, i.e. the according RingList
		boost::property_map<P_Graph, boost::edge_name_t>::type pGraphPath;

		  //! the degree of each node in pGraph
		std::vector<size_t> pGraphDegree;
		  //! the nodes to be removed during the ring perception from pGraph
		std::vector<size_t> toRemove;

	protected:

		  /*!
		   * Initializes the pGraph member with the given graph. If direct loops
		   * are found, these are directly reported and counted.
		   *
		   * @param graph the graph to be encoded in the pGraph member
		   * @param reporter the RingReporter to report all rings to
		   * @return the number of rings found
		   */
		size_t
		initializeP_Graph( const Graph_Interface & graph
							, RingReporter& reporter);

		  /*!
		   * Virtually removes the vertex nextToRemove from the pGraph and
		   * reports all found rings to the given reporter.
		   * The removal is only virtual, since the node remains within the
		   * pGraph but all adjacent edges are removed, i.e. it is not
		   * accessible anymore in later iterations. This is to ease the
		   * maintenance of the temporary data structures.
		   *
		   * @param nextToRemove the index of the vertex to be removed
		   * @param graph the graph searched (needed for ring report)
		   * @param reporter the RingReporter to report all rings to
		   * @param maxRingSize the maximal ring size to consider
		   * @return the number of rings found
		   */
		size_t
		remove_vertex( const size_t nextToRemove
						, const Graph_Interface& graph
						, RingReporter& reporter
						, const size_t maxRingSize);

		struct degree_sort {
			const std::vector<size_t> & degree;
			degree_sort( const std::vector<size_t> & degree_ )
			 : degree(degree_)
			{}
			bool operator()(const size_t& i1, const size_t& i2 ) const {
				return degree[i1] > degree[i2];
			}
		};


		/**
		 * Stores the bond information for loops.
		 */
		class LoopBondReporter : public RingReporter {

		protected:

			BondSet & bondSet;

		public:
			  //! construction
			LoopBondReporter( BondSet & bondSet );
			  //! destruction
			virtual ~LoopBondReporter();

			  /*!
			   * Is called to report a ring.
			   * @param graph the graph that contains the ring
			   * @param ringList the ring to report
			   */
			virtual
			void
			reportRing( const Graph_Interface& graph, const RingList & ringList );

		};

	};

} // namespace sgm

#endif /* SGM_RP_HANSER96_HH_ */
