
#ifndef GGL_CHEM_AP_NSPDK_MODEL_HH_
#define GGL_CHEM_AP_NSPDK_MODEL_HH_

#include <string>
#include <map>
#include <set>

#include <boost/shared_ptr.hpp>


//  // set spirit closure limit if neccessary
//#if !defined(BOOST_SPIRIT_CLOSURE_LIMIT)
//#define BOOST_SPIRIT_CLOSURE_LIMIT 5
//#elif BOOST_SPIRIT_CLOSURE_LIMIT < 5
//#error "GGL_GRAPH_GML_GRAMMAR : BOOST_SPIRIT_CLOSURE_LIMIT too low, has to be at least 5"
//#endif
//
//  // set phoenix limit if neccessary
//#if !defined(PHOENIX_LIMIT)
//#define PHOENIX_LIMIT 5
//#elif PHOENIX_LIMIT < 5
//#error "GGL_GRAPH_GML_GRAMMAR : PHOENIX_LIMIT too low, has to be at least 5"
//#endif
//
//#include <boost/version.hpp>
//#if BOOST_VERSION >= 103800
//#include <boost/spirit/include/qi.hpp>
//#include <boost/spirit/include/phoenix_core.hpp>
//#include <boost/spirit/include/classic.hpp>
//#include <boost/spirit/include/phoenix1.hpp>
//#include <boost/spirit/include/classic_actor.hpp>
//#define NS_BOOSTSPIRIT boost::spirit::classic
//#else
//#include <boost/spirit.hpp>
//#include <boost/spirit/phoenix.hpp>
//#include <boost/spirit/actor.hpp>
//#define NS_BOOSTSPIRIT boost::spirit
//#endif


#include "sgd/svmmodel.h"

#include "nspdk/SVector.h"


namespace ggl {
namespace chem {

	 /*! @brief NSPDK aromaticity model
	  *
	  * A linear SVM model derived using the NSPDK kernel that defines an
	  * aromaticity classifier used by ggl::chem::AP_NSPDK to identify
	  * aromatic rings.
	  *
	  * @author Martin Mann (c) 2011  http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class AP_NSPDK_Model : public sgd::SvmModel {

	public:
		  //! construction
		AP_NSPDK_Model();
		  //! destruction
		~AP_NSPDK_Model();

		  //! the ID of the model
		std::string modelID;

		  //! a description of the model, e.g. its sources etc.
		std::string description;

		  //! whether or not the NSPDK features are ring centered or not
		bool ringCentered;

		  /*!
		   * If ringCentered == true, this parameter defines up to what
		   * distance of ring nodes nodes are viewPoints, ie. a value of 0
		   * restricts the viewpoints to the ring only.
		   */
		size_t nspdk_maxRingDistance;

		  //! the NSPDK distance parameter used to generate the model
		size_t nspdk_maxDistance;

		  //! the NSPDK radius parameter used to generate the model
		size_t nspdk_maxRadius;

		  //! the NSPDK parameter that defines the number of bits that
		  //! define the feature space size used to generate the model
		size_t nspdk_featureBitSize;

		  //! the maximal ring size considered when the model was created
		size_t maxRingSize;

		  //! the prefix added to all node label within all rings,
		  //! used to minimize the uncertain information for prediction
		  //! NOTE: the node label has to be set non-aromatic as well!
		std::string ringNodeLabelPrefix;

		  //! the edge label used within all rings,
		  //! used to minimize the uncertain information for prediction
		std::string ringEdgeLabel;

		  //! the label prefix that is added in front of each node and edge
		  //! label participating in the ring currently under consideration
		std::string ringViewLabelPrefix;

		  //! the weights for each feature
		using sgd::SvmModel::w;

		  //! the weight scale for better precision
		using sgd::SvmModel::wscale;

		  //! the bias that marks the treshold to reach to be classified as a
		  //! positive instance
		using sgd::SvmModel::bias;


	protected : // static members

		typedef boost::shared_ptr<AP_NSPDK_Model> AP_NSPDK_Model_Pointer;
		  //! container of all available predefined models
		static std::map < std::string, AP_NSPDK_Model_Pointer > availableModels;
		  //! list of all available IDs of predefined models
		static std::set < std::string > availableModelIDs;

	public: // static functions

		  /*!
		   * Access to the predefined instances via their identifier.
		   *
		   * @param modelID the identifier of the model of interest
		   * @return the model or NULL if no model is available for the given
		   *         identifier
		   */
		static
		const AP_NSPDK_Model * const
		getInstance( const std::string & modelID );

		  /*!
		   * Enables the access to all registered predefined model IDs that are
		   * handled via getInstance(modelID).
		   * @return the set of all available model IDs
		   */
		static
		const std::set < std::string > &
		getAvailableModels( void );

	protected:  // static model data

		  //! the number of entries covered by the model "OpenBabel:2013"
		static int OpenBabel_2013_wdatasize;
		  //! the feature weight entries of the model "OpenBabel:2013"
		static sgd::SVector::Pair OpenBabel_2013_wdata[];

		  //! the number of entries covered by the model "Marvin:general:2013"
		static int Marvin_general_2013_wdatasize;
		  //! the feature weight entries of the model "Marvin:general:2013"
		static sgd::SVector::Pair Marvin_general_2013_wdata[];


	};


}} // namespace

#endif /* GGL_CHEM_AP_NSPDK_MODEL_HH_ */
