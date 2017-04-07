#ifndef GGL_CHEM_AROMATICITYPERCEPTION_HH_
#define GGL_CHEM_AROMATICITYPERCEPTION_HH_

#include <utility>
#include <set>
#include <stdexcept>

#include "sgm/RingReporter.hh"

#include "ggl/chem/Molecule.hh"

namespace ggl {
namespace chem {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



	 /*! @brief Interface - aromaticity prediction
	  *
	  * Interface to identify aromatic rings in molecules and to relabel it
	  * accordingly.
	  *
	  *  @author Martin Mann (c) 2011 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class AromaticityPerception : public sgm::RingReporter {

	public:

		  //! defines an edge; NOTE: use as ordered pairs, i.e. first <= second
		typedef std::pair< size_t, size_t > Edge;
		  //! container to maintain aromatic edges without duplicates
		typedef std::set< Edge > EdgeSet;


		  /*!
		   * Describes a ring based on a sorted set of edges and the list of
		   * nodes to encode the ring direction.
		   */
		class RingDescriptor {
		public:

			enum AromState { Unknown, Aromatic, NonAromatic };

			  //! edges participating in the ring
			EdgeSet edges;
			  //! nodes participating in the ring in ring order (note: the
			  //! reverse order is not stored but valid as well)
			RingList ring;

			  //! the predicted aromaticity state of this ring
			AromState predState;

			  //! a value >= 0 describing the certainty of the given predState.
			  //! the larger the value the more certain the prediction is.
			double predCertainty;

			  //! empty construction
			RingDescriptor();

			  //! construction based on a ring list encoding
			  //! @ring the ring to be described
			RingDescriptor( const RingList& ring );

		};

	protected:

		  //! the list of edges part of aromatic rings to be identified
		EdgeSet aromaticEdges;

		  //! container that will hold all rings of the current molecule
		  //! -> this list is later pruned to the rings of interest
		std::vector< RingDescriptor* > allRings;

		  /*! Describes the edge information of a node participating in a ring
		   * to enable the correct relabeling of ring edges.
		   *
		   */
		struct AdjacencyData {
		public:
			  //! ID of the node described
			size_t nodeID;
			  //! the remaining valence of the node to be distributed among the
			  //! adjacent ring edges, i.e. the valence of the atom minus the
			  //! bonds valence of all labeled adjacent bonds
			size_t remValence;
			  //! the number of non-labeled adjacent ring edges
			size_t openEdges;
			  //! the number of adjacent labeled aromatic ring edges
			size_t aromaticEdges;
			 /*! Construction
			  * @param nID node ID
			  * @param val remaining valence to be distributed
			  * @param openEdges number of adjacent non-labeled ring edges to relabel
			  * @param aromEdges number of adjacent labeled aromatic ring edges
			  */
			AdjacencyData( 	const size_t nID = 0,
							const size_t val = 0,
							const size_t openEdges = 0,
							const size_t aromEdges = 0 )
			 : nodeID(nID), remValence(val), openEdges(openEdges), aromaticEdges(aromEdges)
			{}


		};

		//! @brief Comparator for edge processing order
		struct AdjacencyComp {
			  /*!
			   * Decides on the order of the edge processing when used to sort
			   * a container of AdjacencyData objects.
			   *
			   * @param e1 the first edge to be checked if to be processed first
			   * @param e2 the second edge that might be processed after e1
			   *
			   * @return whether or not the edge e1 should be processed before
			   *         edge e2
			   */
			bool operator() ( AdjacencyData* e1, AdjacencyData* e2);
		};


		  /*!
		   * Comparison of ring descriptors based on the ring size and if equal
		   * on the pairwise order of its edge elements.
		   */
		class RingSizeLess {
		public:

			  /*!
			   * less comparison of two ring descriptors.
			   * @param x the ring to be checked if smaller
			   * @param y the ring to compare to
			   * @return true if x.size < y.size
			   * 		|| (x.size == y.size && x.edges < y.edges);
			   * 		false otherwise
			   */
			bool operator () ( const RingDescriptor* x
							  , const RingDescriptor* y )
			{
				return x->edges.size() < y->edges.size()
						|| (x->edges.size() == y->edges.size() && x->edges < y->edges );
			}
		};


	public:

		 /*! construction
		  */
		AromaticityPerception();

		 /*! copy construction
		  */
		AromaticityPerception( const AromaticityPerception & toCopy );

		 /*!
		  * destruction
		  */
		virtual ~AromaticityPerception();


		 /*!
		  * Identifies aromatic rings within the molecule and relabels atoms
		  * and bonds accordingly. If rings are wrongly assigned aromatic, they
		  * are relabeled as well.
		  *
		  * @param mol the molecule to correct if needed (done inplace)
		  * @param checkValence whether or not the valence of nodes and bonds
		  *        should be checked. In error case a std::runtime_error
		  *        is thrown.
		  *
		  * @throw std::runtime_error if a canonical relabeling is not possible
		  *          or the valence checks fail
		  */
		virtual
		void
		correctAromaticity( Molecule & mol, const bool checkValence )
			throw (std::runtime_error);

		 /*!
		  * Creates a heap copy of this instance that has to be deleted by the
		  * calling methods later on.
		  * @return a new instance that is a copy of this object
		  */
		virtual
		AromaticityPerception *
		clone() const = 0;

		  /*!
		   * Is called to report a ring identified by a RingPerception instance.
		   *
		   * Here it is used to store each ring within the allRings container.
		   *
		   * @param graph the graph that contains the ring
		   * @param ringList the ring to report
		   */
		virtual
		void
		reportRing(	const sgm::Graph_Interface& graph
					, const RingList & ringList );



	protected:

		  /*!
		   * Resets internal temporary data structures.
		   */
		void
		clearData();

		  /*!
		   * Identifies all rings and stores them within the allRings container.
		   *
		   * This function is abstract and has to be implemented by a
		   * specific subclass. An according RingReporter function is part of
		   * this class.
		   *
		   * @param mol the molecule to check for rings
		   */
		virtual
		void
		findAllRings( const Molecule & mol ) = 0;

		  /*!
		   * Identifies aromatic rings and stores their edges within the
		   * aromaticEdges container.
		   *
		   * This function is abstract and has to be implemented by a
		   * specific subclass.
		   *
		   * @param mol the molecule to check for aromatic rings
		   */
		virtual
		void
		identifyAromaticEdges( const Molecule & mol ) = 0;

		  /*!
		   * Relabels the given molecule according to the aromatic rings
		   * defined by the aromaticEdges container. All other nodes
		   * and edges are set non-aromatic.
		   *
		   * @param mol the molecule to relabel (done inplace)
		   * @param checkValence whether or not the valence of nodes and bonds
		   *        should be checked. In error case a std::runtime_error
		   *        is thrown.
		   *
		   * @throw std::runtime_error if a canonical relabeling is not possible
		   *          or the valence checks fail
		   */
		void
		relabelMolecule( Molecule & mol
						, const bool checkValence ) throw (std::runtime_error);


		  /*!
		   * Performs a pruning of rings (within the "allRings" container)
		   * that are fused versions of smaller rings sharing only one bond.
		   * @param rings the rings to be pruned
		   */
		static
		void
		pruneFusedRings( std::vector< RingDescriptor* > & rings );

		  /*!
		   * Performs a pruning of rings (within the "allRings" container)
		   * that are not composed of single and double bonds only.
		   * @param rings the rings to be pruned
		   * @param mol the molecule to relabel
		   */
		static
		void
		pruneNonSingleDoubleBondRings( std::vector< RingDescriptor* > & rings
							, const Molecule & mol );

	};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} } // namespace

// include implementation
#include "ggl/chem/AromaticityPerception.icc"

#endif /* GGL_CHEM_AROMATICITYPERCEPTION_HH_ */
