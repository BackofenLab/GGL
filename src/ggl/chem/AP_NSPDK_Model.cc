
#include "ggl/chem/AP_NSPDK_Model.hh"

#include <sstream>

#include <boost/algorithm/string.hpp>




namespace ggl {
namespace chem {

	std::map< std::string, AP_NSPDK_Model::AP_NSPDK_Model_Pointer >
		AP_NSPDK_Model::availableModels;

	std::set< std::string >
		AP_NSPDK_Model::availableModelIDs;


	AP_NSPDK_Model::
	AP_NSPDK_Model()
	{
	}

	AP_NSPDK_Model::
	~AP_NSPDK_Model()
	{
	}

	const std::set < std::string > &
	AP_NSPDK_Model::
	getAvailableModels( void )
	{
		  // initialize list if not done yet
		if (availableModelIDs.empty()) {
			  // OpenBabel model
			availableModelIDs.insert("OpenBabel:2013");
			  // Marvin general model
			availableModelIDs.insert("Marvin:general:2013");
		}
		  // access to list
		return availableModelIDs;
	}

///////////////////////////////////////////////////////////////////////////////

	const AP_NSPDK_Model * const
	AP_NSPDK_Model::
	getInstance( const std::string & modelID )
	{

		  // check if the available models are already initialized
		if (availableModels.find(modelID) == availableModels.end()) {

			// do model initialization on demand

			////////////////////////////////////////////////////////////////////
			//////////   "OpenBabel:2013"  ////////////////////////
			////////////////////////////////////////////////////////////////////
			if (modelID == "OpenBabel:2013"){
				  // create new model
				AP_NSPDK_Model* aromaticityModel = new AP_NSPDK_Model() ;
				  // set model identifier
				aromaticityModel->modelID = modelID;
				  // include model definition
				#include "ggl/chem/AP_NSPDK_Model_OpenBabel_2013.icc"
				  // load weight data
				sgd::SVector::Rep r;
				r.npairs = r.mpairs = OpenBabel_2013_wdatasize;
				r.pairs = OpenBabel_2013_wdata;
				r.size = r.pairs[r.npairs-1].i + 1;
				aromaticityModel->w = r;
				 // undo to prevent deletion of static data via ~Rep
				r.npairs = r.mpairs = r.size = 0;
				r.pairs = 0;
				  // store
				availableModels[modelID] = AP_NSPDK_Model_Pointer(aromaticityModel);

			}///////////////////////////////////////////////////////////////////
			else
			////////////////////////////////////////////////////////////////////
			//////////   "Marvin:general:2013"  ////////////////////////
			////////////////////////////////////////////////////////////////////
			if (modelID == "Marvin:general:2013"){
				  // create new model
				AP_NSPDK_Model* aromaticityModel = new AP_NSPDK_Model() ;
				  // set model identifier
				aromaticityModel->modelID = modelID;
				  // include model definition
				#include "ggl/chem/AP_NSPDK_Model_Marvin_general_2013.icc"
				  // load weight data
				sgd::SVector::Rep r;
				r.npairs = r.mpairs = Marvin_general_2013_wdatasize;
				r.pairs = Marvin_general_2013_wdata;
				r.size = r.pairs[r.npairs-1].i + 1;
				aromaticityModel->w = r;
				 // undo to prevent deletion of static data via ~Rep
				r.npairs = r.mpairs = r.size = 0;
				r.pairs = 0;
				  // store
				availableModels[modelID] = AP_NSPDK_Model_Pointer(aromaticityModel);

			}///////////////////////////////////////////////////////////////////
			else {

				  // model not known
				return NULL;
			}
		}
		  // get access to the requested model
		return availableModels.find(modelID)->second.get();
	}

///////////////////////////////////////////////////////////////////////////////

}} // namespaces

