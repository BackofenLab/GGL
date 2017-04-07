
#ifndef GGL_CHEM_AP_NSPDK_HH_
#define GGL_CHEM_AP_NSPDK_HH_

#include <cassert>
#include <vector>
#include <set>
#include <algorithm>
#include <limits>

#include "sgd/svmmodel.h"

#include "nspdk/GraphClass.h"

#include "sgm/Graph_Interface.hh"
#include "sgm/NSPDK_port.hh"

#include "ggl/chem/AromaticityPerception.hh"
#include "ggl/chem/AP_NSPDK_Model.hh"

namespace ggl {
namespace chem {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


	 /*! @brief NSPDK-based aromaticity prediction
	  *
	  * Implementation of an AromaticityPerception based on the NSPDK graph
	  * kernel and linear SVM models derived from chemical databases. The graph
	  * kernel features are derived are all anchored within the according ring
	  * of interest.
	  *
	  *  @author Martin Mann (c) 2011 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class AP_NSPDK
		: public AromaticityPerception
	{
	public :

		  //! ring describing class
		typedef AromaticityPerception::RingDescriptor RingDescriptor;

		  //! Vector that holds ( index, weight ) pairs
		typedef std::vector< std::pair<size_t,double> > RingWeightVec;


	protected:

		  //! defines an edge; NOTE: use as ordered pairs, i.e. first <= second
		typedef AromaticityPerception::Edge AromaticEdge;
		  //! defines an aromatic ring via the set of aromatic edges
		typedef AromaticityPerception::EdgeSet AromaticEdgeSet;
		  //! defines a list of view points for NSPDK feature generation
		typedef sgm::NSPDK_port::ViewPointList ViewPointList;

		  //! the aromatic ring container to fill
		using AromaticityPerception::aromaticEdges;

		  //! container that will hold all rings of the current molecule
		  //! -> this list is later pruned to the rings of interest
		using AromaticityPerception::allRings;

		  //! the NSPDK graph instance used for feature generation and testing
		nspdk::GraphClass nspdkGraph;

		  //! the SVM model used to predict the aromaticity of a ring
		const AP_NSPDK_Model* aromaticityModel;


	public:

		  /*!
		   * Construction
		   * @param modelID the identifier of the aromaticity model to be used
		   *        (check ggl::chem::AP_NSPDK_Model)
		   */
		AP_NSPDK( const std::string & modelID );

		  /*!
		   * Construction
		   * @param model the aromaticity model to be used
		   */
		AP_NSPDK( const AP_NSPDK_Model * const model );

		  /*!
		   * Copy construction
		   * @param toCopy the object to make this a copy of
		   */
		AP_NSPDK( const AP_NSPDK & toCopy );

		  /*!
		   * destruction
		   */
		virtual ~AP_NSPDK();


		  /*!
		   * Identifies all rings up to the maximally predictable size of the
		   * used NSPDK model and stores them within the allRings container.
		   *
		   * @param mol the molecule to check for rings
		   */
		virtual
		void
		findAllRings( const Molecule & mol );

		  /*!
		   * Identifies aromatic rings and stores their edges within the
		   * aromaticEdges container.
		   *
		   * To this end, each ring is classified using an linear NSPDK SVM
		   * model of aromaticity.
		   *
		   * @param mol the molecule to check for aromatic rings
		   */
		void
		identifyAromaticEdges( const Molecule & mol );


		  /*!
		   * Access to the aromaticity model used.
		   */
		const AP_NSPDK_Model* const
		getModel( void ) const;

		 /*!
		  * Creates a heap copy of this instance that has to be deleted by the
		  * calling methods later on.
		  * @return a new instance that is a copy of this object
		  */
		virtual
		AP_NSPDK *
		clone() const;


		  /*!
		   * Identifies all rings and returns for each a RingList representation
		   * as well as the according NSPDK features.
		   *
		   * Note: if the used model is ringCentered, only ring nodes are
		   * viewpoints (ring.ring.begin()..(ring.ring.end()-1)))
		   * and thus returned; otherwise a weight for each node is
		   * returned.
		   *
		   * @param mol the molecule to check for rings
		   * @return a ring representation and the features of each ring
		   */
		virtual
		std::vector< std::pair< RingList, nspdk::SVector > >
		getFeatures( const Molecule & mol );

		  /*!
		   * Identifies all rings and returns for each a RingDescriptor
		   * as well as the normalized feature weight for each viewpoint node.
		   *
		   * Note: if the used model is ringCentered, only ring nodes are
		   * viewpoints (ring.ring.begin()..(ring.ring.end()-1)))
		   * and thus returned; otherwise a weight for each node is
		   * returned.
		   *
		   * @param mol the molecule to check for rings
		   * @return a ring representation and the weights for each ring node
		   */
		virtual
		std::vector< std::pair< RingDescriptor, RingWeightVec > >
		getNodeWeights( const Molecule & mol );

	protected:  // member functions

		  /*!
		   * Computes the NSPDK feature vector for a given ring within a
		   * molecule. Note, the molecule graph has to show already the
		   * according ring relabeling.
		   *
		   * Note: if the used model is ringCentered, only ring nodes are
		   * viewpoints (ring.ring.begin()..(ring.ring.end()-1)))
		   * and thus returned; otherwise a weight for each node is
		   * returned.
		   *
		   * @param molGraph the molecule graph containing the ring
		   * @param ring the ring to be described
		   * @param aromaticityModel the model to use
		   */
		virtual
		nspdk::SVector
		getFeaturesOfRing( nspdk::GraphClass &molGraph
							, const RingDescriptor & ring
							, const AP_NSPDK_Model& aromaticityModel);

		  /*!
		   * Computes the normalized NSPDK feature weights for each viewpoint
		   * node in the graph for a ring within a molecule.
		   * Note, the molecule graph has to show already the
		   * according ring relabeling.
		   *
		   * Note further: if the model is ringCentered, only ring nodes are
		   * viewpoints (ring.ring.begin()..(ring.ring.end()-1)))
		   * and thus returned; otherwise a weight for each node is
		   * returned.
		   *
		   * @param molGraph the molecule graph containing the ring
		   * @param ring the ring to be described
		   * @param aromaticityModel the model to use
		   * @return the feature vector for each viewpoint (see above)
		   */
		virtual
		RingWeightVec
		getWeightsOfRing( nspdk::GraphClass &molGraph
							, const RingDescriptor & ring
							, const AP_NSPDK_Model& aromaticityModel );

		  /*!
		   * Does an aromaticity prediction for the given ring.
		   *
		   * @param ring the ring to predict
		   * @return true if the ring is aromatic; false otherwise
		   */
		virtual
		bool
		isAromaticRing( const RingDescriptor & ring );

		  /*!
		   * Decides on an aromaticity prediction value if it means aromatic
		   * or not.
		   *
		   * @param predictionValue the prediction value to decide on
		   * @return true if the ring is aromatic; false otherwise
		   */
		virtual
		bool
		isAromaticRing( const double predictionValue );

		  /*!
		   * Computes the classification prediction value for this ring.
		   * This value can be given to "isAromaticRing()"
		   * in order to decide if it is aromatic or not.
		   *
		   * @param ring the ring to predict
		   * @return the prediction value for this ring
		   */
		virtual
		double
		getPredictionValue( const RingDescriptor & ring );

	protected:  // member classes


		  /*!
		   * Generates the set of view point nodes for a given ring in the
		   * molecule graph that anchor the feature generation for the NSPDK
		   * graph kernel. The view point set depends on the aromaticity model
		   * setup, ie. if its covering all nodes or only the ring nodes and
		   * their surrounding.
		   *
		   * @param molGraph the graph of interest
		   * @param ring the ring of interest
		   * @param model the model defining what view points to choose
		   */
		virtual
		ViewPointList
		getViewPoints( nspdk::GraphClass & molGraph
						, const RingDescriptor & ring
						, const AP_NSPDK_Model & model ) const;


		  /*!
		   * Initializes the nspdkGraph and relabels all nodes and edges
		   * of the nspdkGraph involved within rings.
		   *
		   * NOTE: rings that contain edges with valence > 2 or < 1 are
		   * ignored!
		   *
		   * @param mol the molecule graph to be used for initialization
		   * @param edgeLabel the edge label to set for all rings
		   * @param ringFlag the flag that is added to the beginning of
		   *         each node label of nodes involved in a ring
		   */
		virtual
		void
		initializeGraph( const Molecule mol
						, const std::string edgeLabel
						, const std::string ringFlag );


		//////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////
	};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}} // namespace

#include "ggl/chem/AP_NSPDK.icc"

#endif /* GGL_CHEM_AP_NSPDK_HH_ */
