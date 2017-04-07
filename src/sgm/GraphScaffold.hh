#ifndef SGM_GRAPHSCAFFOLD_HH_
#define SGM_GRAPHSCAFFOLD_HH_


#include "sgm/Graph_Interface.hh"
#include "sgm/RP_Hanser96.hh"

namespace sgm {


	/*!
	 * Computes the scaffold annotation of each node, ie. whether it is part of
	 * a ring, a ring linker, or a dangling end (side chain).
	 *
	 * See
	 *   The Properties of Known Drugs. 1. Molecular Frameworks
	 *   Guy W. Bemis and Mark A. Murcko
	 *   J. Med. Chem. 1996, 39, 2887-2893
	 *
	 * @author Martin Mann (c) 2013 http://www.bioinf.uni-freiburg.de/~mmann/
	 *
	 */
	class GraphScaffold : public RP_Hanser96 {

	public:

		  //! possible scaffold types of graph nodes
		enum ScaffoldType {
			GST_UNKNOWN,	//!< nodes of this type are not processed yet
			GST_RING, 	//!< nodes of this type are part of a ring
			GST_LINKER,	//!< nodes of this type are part of a ring linker
			GST_DANGLING	//!< nodes of this type are part of a dangling end / side chain
		};

		  //! assignment of scaffold types for each node index
		typedef std::vector< ScaffoldType > ScaffoldAnnotation;

	protected:

		  //! the P-graph to be compressed during the ring perception
		using RP_Hanser96::pGraph;

		  //! the index access for pGraph
		using RP_Hanser96::pGraphIndex;
		  //! the edge label access for pGraph, i.e. the according RingList
		using RP_Hanser96::pGraphPath;

		  //! the degree of each node in pGraph
		using RP_Hanser96::pGraphDegree;
		  //! the nodes to be removed during the ring perception from pGraph
		using RP_Hanser96::toRemove;


	public:

		  //! construction
		GraphScaffold();

		  //! destruction
		virtual ~GraphScaffold();

		  /*!
		   * Identifies the scaffold type for each node in the graph and returns
		   * the according scaffold annotation.
		   *
		   * @param graph the graph to annotate
		   * @param maxRingSize the maximal size of rings to consider
		   * @return the scaffold type for each node within the graph, ie. the
		   * container has as many entries as the graph nodes.
		   */
		ScaffoldAnnotation
		getScaffoldAnnotation( const Graph_Interface & graph
				, const size_t maxRingSize = std::numeric_limits<size_t>::max() );

		  /*!
		   * Identifies the scaffold type for each node in the graph and returns
		   * the according scaffold annotation.
		   *
		   * @param graph the graph to annotate
		   * @param reporter the ring reporter to forward all rings to
		   * @param maxRingSize the maximal size of rings to consider
		   * @return the scaffold type for each node within the graph, ie. the
		   * container has as many entries as the graph nodes.
		   */
		ScaffoldAnnotation
		getScaffoldAnnotation( const Graph_Interface & graph
				, RingReporter & reporter
				, const size_t maxRingSize = std::numeric_limits<size_t>::max());

	protected:




		  /*!
		   * Identifies the scaffold type for each node in the graph and returns
		   * the according scaffold annotation.
		   *
		   * @param graph the graph to annotate
		   * @param reporter the ring reporter to forward all rings to (can be
		   * NULL such that no reporting is done)
		   * @param maxRingSize the maximal size of rings to consider
		   * @return the scaffold type for each node within the graph, ie. the
		   * container has as many entries as the graph nodes.
		   */
		ScaffoldAnnotation
		getScaffoldAnnotation( const Graph_Interface & graph
								, RingReporter * reporter
								, const size_t maxRingSize );


		 /*!
		  * Ring reporter that stores ring annotations and forwards the
		  * reporting to another ring reporter if given.
		  */
		class RR_Annotation : public RingReporter {

		protected:

			  //! the annotation to update with ring information
			ScaffoldAnnotation & annotation;
			  //! the next reporter to forward the ring information to if not NULL
			RingReporter * nextReporter;

		public:

			 /*!
			  * Constructions
			  * @param annotation the annotation to update with ring information
			  * @param nextReporter the next reporter to forward the ring
			  * information to if not NULL
			  */
			RR_Annotation( ScaffoldAnnotation & annotation
						, RingReporter * nextReporter = NULL)
			 : annotation(annotation)
				, nextReporter(nextReporter)
			{}

			virtual
			~RR_Annotation()
			{}
			  /*!
			   * Nodes of the reported ring are accordingly annotated and the
			   * ring information is forwarded to the next ring reporter
			   * @param graph the graph that contains the ring
			   * @param ringList the ring to report
			   */
			virtual
			void
			reportRing( const Graph_Interface& graph, const RingList & ringList );

		};

	};

} // namespace sgm

#endif /* SGM_GRAPHSCAFFOLD_HH_ */
