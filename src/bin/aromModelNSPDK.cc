/*
 * This binary enables the training and testing of aromaticity models using
 * the NSPDK kernel within the SGD framework.
 *
 *      Author: mmann
 */



#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <exception>
#include <iterator>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include <map>

#include <boost/algorithm/string.hpp>

#include "biu/OptionParser.hh"
#include "biu/Timer.hh"

#include <sgm/Graph_boost.hh>

#include <ggl/Graph.hh>
#include <ggl/chem/Molecule.hh>
#include <ggl/chem/AP_NSPDK.hh>
#include <ggl/chem/SMILESparser.hh>
#include <ggl/chem/SMILESwriter.hh>
#include <ggl/Graph_GML_writer.hh>

#include <sgd/svmsgd.h>

#include "version.hh"

	using namespace ggl;
	using namespace ggl::chem;


//////////////////////////////////////////////////////////////////////////
#ifndef ARGEXCEPTION_
#define ARGEXCEPTION_

	  /*! Exception class for exeptions thrown during argument and input parsing.
	   */
	class ArgException : public std::exception {
	public:
		  //! the error message
		std::string errorMsg;
		ArgException( std::string errorMsg_ ) : errorMsg(errorMsg_) {}
		virtual ~ArgException() throw (){}
		virtual const char* what() const throw() {
			return errorMsg.c_str();
		}
	};

#endif
//////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

bool
SVector_less(const sgd::SVector& v, const sgd::SVector & vCompare);

/**
 * Dummy pointer class that represents an index within a vector and enables
 * vector element comparison/sorting without vector alteration.
 */
class FeatureInfo {

	size_t posInVector;
	sgd::xvec_t * features;

public:

	FeatureInfo()
	: posInVector(0), features(NULL)
	{}

	FeatureInfo( size_t posInVector
				, sgd::xvec_t * features )
	: posInVector(posInVector)
		, features(features)
	{}

	size_t getPosition() const {
		return posInVector;
	}

	bool
	operator < ( const FeatureInfo & toComp ) const {
		assert( toComp.features == this->features );
		return SVector_less( this->features->at(posInVector), toComp.features->at(toComp.posInVector) );
	}
};

/////////////////////////////////////////////////////////////////////////////

/**
 * IN CLUSTERING MODE : Collects information for each ring cluster.
 */
class ClusterInfo {
public:
	  //! position of the representative within the feature vector
	size_t posInVector;
	  //! number of features for this cluster
	size_t clusterSize;
	  //! number of molecules with this feature
	size_t molNumber;

	ClusterInfo()
	  : posInVector(0), clusterSize(0), molNumber(0)
	{}

	ClusterInfo( const size_t posInVector, const size_t clusterSize, const size_t molNumber)
	  : posInVector(posInVector), clusterSize(clusterSize), molNumber(molNumber)
	{}

};

/////////////////////////////////////////////////////////////////////////////

void
annotateRing( Molecule & mol, const sgm::RingReporter::RingList & ring, bool addAnnotation );

/////////////////////////////////////////////////////////////////////////////

void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );

//////////////////////////////////////////////////////////////////////////

void getFeatures( AP_NSPDK & ap
				, const Molecule & mol
				, sgd::xvec_t& trainFeatures
				, sgd::yvec_t& trainTarget
				, std::vector< sgm::RingReporter::RingList > * trainRing = NULL
				);

//////////////////////////////////////////////////////////////////////////

void getFeatures(	AP_NSPDK & ap
					, const std::string & molSMILES
					, sgd::xvec_t& trainFeatures
					, sgd::yvec_t& trainTarget );

/////////////////////////////////////////////////////////////////////////////

typedef std::vector< std::string > SMILES_container;
typedef std::back_insert_iterator< SMILES_container > SMILES_inserter;

void parseSMILES(	std::istream & in
					, SMILES_inserter & toFill
					, const size_t linesRead ) throw(std::exception);

/////////////////////////////////////////////////////////////////////////////

void
readModel( std::istream & in, AP_NSPDK_Model & model );

/////////////////////////////////////////////////////////////////////////////

void
printModel( std::ostream & out, const AP_NSPDK_Model & model );

/////////////////////////////////////////////////////////////////////////////

void
printWeightedGML( std::ostream & out, AP_NSPDK & apNSPDK, const Molecule & mol, const bool multilineGML );

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

enum Mode { MODE_TRAIN, MODE_TEST, MODE_APPLY, MODE_CLUST, MODE_XVAL, MODE_UNDEF };

int main( int argc, char** argv ) {

	//////////////////////////////////////////////////////////////
	// variables
	//////////////////////////////////////////////////////////////

	std::istream* inSMILES = NULL;
	std::ifstream* inSMILESFile = NULL;
	size_t inSMILESLine = 0;

	std::istream* inModel = NULL;
	std::ifstream* inModelFile = NULL;

	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;

	std::ostream* log = &std::cout;
	std::ofstream* logFile = NULL;

	std::ostream* prediction = NULL;
	std::ofstream* predictionFile = NULL;

	SMILES_container SMILES;
	Mode mode = MODE_UNDEF;

	double lambda = 0.0;
	size_t epochs = 0;

	std::string modelID = "";
	size_t maxRingSize = 0;
	size_t nspdkMaxRingDistance = 0;
	size_t nspdkMaxDistance = 0;
	size_t nspdkMaxRadius = 0;
	size_t nspdkFeatureBitSize = 0;
	size_t randomClusterLists = 0;

	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////

	biu::OptionMap allowedArgs;	//< map of allowed arguments
	std::string infoText;		//< info string of the program
	initAllowedArguments(allowedArgs,infoText);	// init

		// parse programm arguments
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
													argv, infoText);
		// check arguments parseable and all mandatory arguments given
		// + help output if needed
	if (opts.argExist("help")) {
		opts.coutUsage();
		return 0;
	}
	if (opts.argExist("version")) {
		giveVersion();
		return 0;
	}
	if (!opts.noErrors()) {
		return -1;
	}

	int exitValue = 0;

	try {

		  // set SMILES input stream
		if (boost::iequals(opts.getStrVal("smiles"),"STDIN")) {
			  // read from STDIN
			inSMILES = &std::cin;
		} else {
			  // read from file
			inSMILESFile = new std::ifstream( opts.getStrVal("smiles").c_str()
										, std::ifstream::in );
			if (!inSMILESFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open SMILES input file '" <<opts.getStrVal("smiles") <<"'";
				throw ArgException(oss.str());
			}
			inSMILES = inSMILESFile;
		}
		  // set output stream
		if (opts.getStrVal("out").size() == 0) {
			throw ArgException("no output file given");
		} else if ( !boost::iequals(opts.getStrVal("out"),"STDOUT")) {
			outFile = new std::ofstream(	opts.getStrVal("out").c_str()
											, std::ofstream::out );
			if (!outFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open output file '" <<opts.getStrVal("out") <<"'";
				throw ArgException(oss.str());
			}
			out = outFile;
		}
		  // set logfile stream
		if (opts.getStrVal("log").size() == 0) {
			throw ArgException("no logfile file given");
		} else if ( boost::iequals(opts.getStrVal("log"),"STDOUT")) {
			log = &std::cout;
		} else if ( boost::iequals(opts.getStrVal("log"),"STDERR")) {
			log = &std::cerr;
		} else {
			logFile = new std::ofstream(	opts.getStrVal("log").c_str()
											, std::ofstream::out );
			if (!logFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open logfile file '" <<opts.getStrVal("log") <<"'";
				throw ArgException(oss.str());
			}
			log = logFile;
		}
		  // set prediction stream if needed
		if (opts.argExist("prediction") && opts.getStrVal("prediction").size() >= 0) {
			if ( boost::iequals(opts.getStrVal("prediction"),"STDOUT")) {
				prediction = &std::cout;
			} else if ( boost::iequals(opts.getStrVal("prediction"),"STDERR")) {
				prediction = &std::cerr;
			} else {
				predictionFile = new std::ofstream(	opts.getStrVal("prediction").c_str()
												, std::ofstream::out );
				if (!predictionFile->is_open()) {
					std::ostringstream oss;
					oss	<<"cannot open prediction file '" <<opts.getStrVal("prediction") <<"'";
					throw ArgException(oss.str());
				}
				prediction = predictionFile;
			}
		}
		  // set mode
		if (boost::iequals(opts.getStrVal("mode"), "train")) {
			mode = MODE_TRAIN;
		} else if (boost::iequals(opts.getStrVal("mode"), "test")) {
			mode = MODE_TEST;
		} else if (boost::iequals(opts.getStrVal("mode"), "apply")) {
			mode = MODE_APPLY;
		} else if (boost::iequals(opts.getStrVal("mode"), "cluster")) {
			mode = MODE_CLUST;
		} else if (boost::iequals(opts.getStrVal("mode"), "xval")) {
			mode = MODE_XVAL;
		} else {
			std::ostringstream oss;
			oss	<<"mode '" <<opts.getStrVal("mode") <<"' is unknown";
			throw ArgException(oss.str());
		}

		  // set SGD lambda
		lambda = opts.getDoubleVal("lambda");
		  // set number of epochs
		if (opts.getIntVal("epochs") <= 0) {
			throw ArgException("number of epochs has to be >= 1");
		} else {
			epochs = (size_t)opts.getIntVal("epochs");
		}

		  // set model ID
		modelID = opts.getStrVal("modelID");
		  // set max ring size
		if (opts.getIntVal("maxRingSize") <= 2) {
			throw ArgException("maximal ring size has to be >= 3");
		} else {
			maxRingSize = (size_t)opts.getIntVal("maxRingSize");
		}
		  // set the NSPDK maximal distance
		if (opts.getIntVal("nspdkMaxDistance") < 0) {
			throw ArgException("NSPDK maximal distance has to be >= 0");
		} else {
			nspdkMaxDistance = (size_t)opts.getIntVal("nspdkMaxDistance");
		}
		  // set the NSPDK maximal radius
		if (opts.getIntVal("nspdkMaxRadius") < 0) {
			throw ArgException("NSPDK maximal radius has to be >= 0");
		} else {
			nspdkMaxRadius = (size_t)opts.getIntVal("nspdkMaxRadius");
		}
		  // set the NSPDK maximal radius
		if (opts.argExist("nspdkMaxRingDistance")) {
			if (!opts.argExist("ringCentered")) {
				throw ArgException("Giving 'nspdkMaxRingDistance' without 'ringCentered' makes no sense..");
			}
			if (opts.getIntVal("nspdkMaxRingDistance") < 0) {
				throw ArgException("The maximal distance from the ring for feature generation has to be >= 0");
			} else {
				nspdkMaxRingDistance = (size_t)opts.getIntVal("nspdkMaxRingDistance");
			}
		}
		  // set the NSPDK feature bit size
		if (opts.getIntVal("nspdkFeatureBitSize") < 8) {
			throw ArgException("NSPDK maximal radius has to be >= 8");
		} else {
			nspdkFeatureBitSize = (size_t)opts.getIntVal("nspdkFeatureBitSize");
		}


		  // set number of random cluster lists
		if (opts.getIntVal("outClustLists") < 0) {
			throw ArgException("The number of random cluster sampling lists has to be at least 1.");
		} else {
			randomClusterLists = (size_t)opts.getIntVal("outClustLists");
		}

		/////////////////////////////////////////////////////////////
		// parse all SMILES
		/////////////////////////////////////////////////////////////

		  // parse SMILES
		{
			SMILES_inserter inserter(SMILES);
			inSMILESLine = 0;
			parseSMILES( *inSMILES, inserter, inSMILESLine );
		}

		switch (mode) {

		case MODE_TRAIN : {

			/////////////////////////////////////////////////////////////
			// training mode
			/////////////////////////////////////////////////////////////

			  // setup new model to be filled
			AP_NSPDK_Model newModel;

			newModel.ringCentered = opts.argExist("ringCentered");
			newModel.modelID = modelID;
			newModel.maxRingSize = maxRingSize;
			newModel.nspdk_maxRingDistance = nspdkMaxRingDistance;
			newModel.nspdk_maxDistance = nspdkMaxDistance;
			newModel.nspdk_maxRadius = nspdkMaxRadius;
			newModel.nspdk_featureBitSize = nspdkFeatureBitSize;

			newModel.ringEdgeLabel = "*~";
			newModel.ringNodeLabelPrefix = "*";
			newModel.ringViewLabelPrefix = "?";

			  // create AP_NSPDK instance using the model
			AP_NSPDK apNSPDK( & newModel );
			  // create features
			sgd::xvec_t trainFeatures;
			sgd::yvec_t trainTarget;
			biu::Timer timer;
			timer.start();
			for (size_t i=0; i<SMILES.size(); ++i ) {
				getFeatures( apNSPDK, SMILES.at(i), trainFeatures, trainTarget );
			}
			*log <<" time generating    = " <<(size_t)timer.stop() <<" ms" <<std::endl;

			  // create SGD instance
			sgd::SvmSgd trainer( lambda );
			  // run training
			timer.start();
			trainer.train( trainFeatures, trainTarget, 0, trainFeatures.size(), epochs );
			*log <<" time training      = " <<(size_t)timer.stop() <<" ms" <<std::endl;

			  // store resulting model
			newModel.bias = trainer.getModel().bias;
			newModel.wscale = trainer.getModel().wscale;
			newModel.w = trainer.getModel().w;

			  // output model
			printModel( *out, newModel );

			  // print statistics
			size_t countArom = 0;
			for (sgd::yvec_t::const_iterator it = trainTarget.begin();
					it != trainTarget.end(); ++it )
			{
				if (*it > 0)
					countArom++;
			}
			*log
				<<" molecules          = " <<SMILES.size()
				<<"\n overall rings      = " <<trainTarget.size()
				<<"\n aromatic rings     = " <<countArom
				<<" = " <<(double(countArom)/double(trainTarget.size())*100.0) <<" %"
				<<"\n non-aromatic rings = " <<(trainTarget.size()-countArom)
				<<" = " <<(double(trainTarget.size()-countArom)/double(trainTarget.size())*100.0) <<" %"
				<<"\n" <<std::endl;

		}; break;
		case MODE_XVAL : {

			/////////////////////////////////////////////////////////////
			// cross validation mode
			/////////////////////////////////////////////////////////////

			  // set number of chunks to randomly split the data into
			if (opts.getIntVal("xval") < 2) {
				throw ArgException("number of epochs has to be >= 2");
			}
			const size_t xval = (size_t)opts.getIntVal("xval");

			  // setup new model to be filled
			AP_NSPDK_Model newModel;

			newModel.ringCentered = opts.argExist("ringCentered");
			newModel.modelID = modelID;
			newModel.maxRingSize = maxRingSize;
			newModel.nspdk_maxRingDistance = nspdkMaxRingDistance;
			newModel.nspdk_maxDistance = nspdkMaxDistance;
			newModel.nspdk_maxRadius = nspdkMaxRadius;
			newModel.nspdk_featureBitSize = nspdkFeatureBitSize;

			newModel.ringEdgeLabel = "*~";
			newModel.ringNodeLabelPrefix = "*";
			newModel.ringViewLabelPrefix = "?";

			  // create AP_NSPDK instance using the model
			AP_NSPDK apNSPDK( & newModel );
			  // create features
			sgd::xvec_t trainFeatures;
			sgd::yvec_t trainTarget;
			biu::Timer timer;
			timer.start();
			for (size_t i=0; i<SMILES.size(); ++i ) {
//				std::cerr <<" generating for " <<SMILES.at(i) <<"\n";
				getFeatures( apNSPDK, SMILES.at(i), trainFeatures, trainTarget );
			}
			*log <<" time generating    = " <<(size_t)timer.stop() <<" ms" <<std::endl;
			*log <<" ring number        = " <<trainFeatures.size() <<std::endl;

			  // create randomized chunks
			std::vector<size_t> randomTrainIds(trainFeatures.size(),0);
			for (size_t i=0; i<randomTrainIds.size(); ++i) {
				randomTrainIds[i] = i;
			}
			 // initialize random seed:
			std::srand( unsigned (std::time(NULL)) );
			 // shuffle indices
			std::random_shuffle(randomTrainIds.begin(),randomTrainIds.end());
			 // get chunk size rounded to ceiling
			const size_t chunkSize = (size_t)std::ceil((double)randomTrainIds.size() / (double)xval);

			 // generate test annotation for each feature vector via cross-validation
			sgd::yvec_t testValue(trainTarget.size(),-17.0);
			  // create SGD instance
			sgd::SvmSgd trainer( lambda );
			for (size_t chunk = 0; chunk < xval; ++chunk) {

				  // GENERATE TRAIN AND TEST DATA

				size_t chunkStart = chunk*chunkSize;
				sgd::xvec_t curTestFeatures( chunkSize - ( chunk+1 == xval ? (chunkStart+chunkSize) - trainFeatures.size() : 0) );
				  // collect train chunks
				sgd::xvec_t curTrainFeatures( trainFeatures.size() - curTestFeatures.size() );
				sgd::yvec_t curTrainTargets( trainFeatures.size() - curTestFeatures.size() );
				  // collect leading chunks
				size_t curPos = 0;
				for (size_t c=0; c<chunk; ++c) {
					chunkStart = c*chunkSize;
					for (size_t i=0; i<chunkSize; ++i) {
						curTrainFeatures[curPos] = trainFeatures.at(randomTrainIds.at(chunkStart+i));
						curTrainTargets[curPos]  = trainTarget.at(randomTrainIds.at(chunkStart+i));
						++curPos;
					}
				}
				  // collect trailing chunks
				for (size_t c=chunk+1; c<xval; ++c) {
					chunkStart = c*chunkSize;
					for (size_t i=0; i<chunkSize && (chunkStart+i) < trainFeatures.size(); ++i) {
						curTrainFeatures[curPos] = trainFeatures.at(randomTrainIds.at(chunkStart+i));
						curTrainTargets[curPos]  = trainTarget.at(randomTrainIds.at(chunkStart+i));
						++curPos;
					}
				}

				  // GENERATE MODEL

				  // run training
				timer.start();
				trainer.train( curTrainFeatures, curTrainTargets, 0, curTrainFeatures.size(), epochs );
				*log <<chunk <<" time training      = " <<(size_t)timer.stop() <<" ms" <<std::endl;

				  // store resulting model
				newModel.bias = trainer.getModel().bias;
				newModel.wscale = trainer.getModel().wscale;
				newModel.w = trainer.getModel().w;

				  // RUN AND STORE PREDICTION ON TEST
				chunkStart = chunk*chunkSize;
				  // collect test chunk
				double timeTesting = 0;
				for (size_t i=0; i<chunkSize && (chunkStart+i) < trainFeatures.size(); ++i) {
					  // test instance
					timer.start();
					testValue[randomTrainIds.at(chunkStart+i)] = newModel.predictValue( trainFeatures.at(randomTrainIds.at(chunkStart+i)) );
					timeTesting += timer.stop();
				}
				*log <<chunk <<" time testing       = " <<(size_t)timer.stop() <<" ms" <<std::endl;

			}

			 // GATHER STATISTICS

			size_t truePositive  = 0; // correctly aromatic
			size_t trueNegative  = 0; // correctly non-aromatic
			size_t falsePositive = 0; // incorrectly aromatic (is non-aromatic)
			size_t falseNegative = 0; // incorrectly non-aromatic (is aromatic)

			  // compare each target and test value
			  // and print target-test value pairs to out
			*out <<"targetVal\tpredVal\n";
			for (size_t i=0; i<testValue.size(); ++i) {
				*out <<trainTarget.at(i) <<"\t" <<testValue.at(i) <<"\n";
				  // correct prediction
				if (trainTarget.at(i)*testValue.at(i) > 0) {
					if (trainTarget.at(i) > 0) {
						++truePositive;
					} else {
						++trueNegative;
					}
				} // incorrect prediction
				else {
					if (trainTarget.at(i) > 0) {
						++falseNegative;
					} else {
						++falsePositive;
					}
				}
			}

			  // print statistics
			size_t countArom = (truePositive+falseNegative);
			*log
				<<" molecules          = " <<SMILES.size()
				<<"\n overall rings      = " <<trainTarget.size()
				<<"\n aromatic rings     = " <<countArom
				<<" = " <<(double(countArom)/double(trainTarget.size())*100.0) <<" %"
				<<"\n non-aromatic rings = " <<(trainTarget.size()-countArom)
				<<" = " <<(double(trainTarget.size()-countArom)/double(trainTarget.size())*100.0) <<" %"
				<<"\n" <<std::endl
				<<"\n TP (arom)            = " <<truePositive
				<<"\n TN (non-arom)        = " <<trueNegative
				<<"\n FN (arom)            = " <<falseNegative
				<<"\n FP (non-arom)        = " <<falsePositive
				<<"\n" <<std::endl
				;


		}; break;
		case MODE_CLUST : {

			/////////////////////////////////////////////////////////////
			// clustering mode
			/////////////////////////////////////////////////////////////


			 // generate Molecule representation of all SMILES
			std::vector< Molecule > molecules(SMILES.size());
			for (size_t i=0; i<SMILES.size(); ++i) {

				  // get molecule graph
				std::pair< Molecule, int > parse = SMILESparser::parseSMILES( SMILES.at(i) );
				if (parse.second >= 0) {
					throw ArgException("SMILES '"+SMILES.at(i)+"' cannot be parsed.");
				}
				  // copy to container
				MoleculeUtil::copy( parse.first, molecules[i] );
				  // fill protons
				MoleculeUtil::fillProtons( molecules[i] );

			}


			  // setup new model to be filled
			AP_NSPDK_Model newModel;

			newModel.ringCentered = opts.argExist("ringCentered");
			newModel.modelID = modelID;
			newModel.maxRingSize = maxRingSize;
			newModel.nspdk_maxRingDistance = nspdkMaxRingDistance;
			newModel.nspdk_maxDistance = nspdkMaxDistance;
			newModel.nspdk_maxRadius = nspdkMaxRadius;
			newModel.nspdk_featureBitSize = nspdkFeatureBitSize;

			newModel.ringEdgeLabel = "*~";
			newModel.ringNodeLabelPrefix = "*";
			newModel.ringViewLabelPrefix = "?";

			  // create AP_NSPDK instance using the model
			AP_NSPDK apNSPDK( & newModel );
			  // create features
			sgd::xvec_t trainFeatures;
			sgd::yvec_t trainTarget;
			  // the ring information for each feature vector
			std::vector< sgm::RingReporter::RingList > trainRing;
			typedef std::vector< size_t > SizetVec;
			SizetVec trainSMILES;
			biu::Timer timer;
			timer.start();
			for (size_t i=0; i<molecules.size(); ++i ) {
				getFeatures( apNSPDK, molecules.at(i), trainFeatures, trainTarget, &trainRing );
				// store SMILES/molecule ID for each added feature
				trainSMILES.resize(trainFeatures.size(),i);
			}
			*log <<" time generating    = " <<(size_t)timer.stop() <<" ms" <<std::endl;
			timer.start();

			std::vector< FeatureInfo > featureInfo( trainFeatures.size() );
			for (size_t i=0; i<featureInfo.size(); ++i) {
				featureInfo[i] = FeatureInfo(i, &trainFeatures);
			}

			  // sort all features but preserve ID information for the features
			std::sort( featureInfo.begin(), featureInfo.end());

			*log <<" time sorting       = " <<(size_t)timer.stop() <<" ms" <<std::endl;

			  // gather histogram data for cluster size
			  // and output at most two random samples per cluster
			typedef std::map<size_t,size_t> HIST;
			HIST clusterSizeHistogram;
			HIST clusterMolHistogram;
			size_t cStart = 0, cNext=0, cID = 0;


//			std::cerr <<"\n#############  ITERATE  ########## ";

			std::vector< std::vector<ClusterInfo> > selectedRings(randomClusterLists);
			/* initialize random seed: */
			std::srand ( unsigned (std::time(NULL)) );
			std::set<size_t> molIDs;
			while( cNext<featureInfo.size() ) {
				molIDs.clear();
				  // insert index of first cluster molecule
				molIDs.insert(trainSMILES.at(featureInfo.at(cNext).getPosition()));
				  // find end of current cluster
				++cNext;
				while(cNext<featureInfo.size()
						&& ! (featureInfo.at(cStart) < featureInfo.at(cNext)))
				{
					  // insert index of next cluster molecule
					molIDs.insert(trainSMILES.at(featureInfo.at(cNext).getPosition()));
					++cNext;
				}

				// update size histogram information
				size_t clusterSize = cNext-cStart;
				if (clusterSizeHistogram.find(clusterSize) == clusterSizeHistogram.end()) {
					clusterSizeHistogram[clusterSize] = 0;
				}
				clusterSizeHistogram[clusterSize]++;

				  // update molecule number histogram information
				size_t molNumber = molIDs.size();
				if (clusterMolHistogram.find(molNumber) == clusterMolHistogram.end()) {
					clusterMolHistogram[molNumber] = 0;
				}
				clusterMolHistogram[molNumber]++;

				  // pick random element of this cluster
				for (size_t i=0; i<randomClusterLists; ++i) {
					if (clusterSize == 1) {
						selectedRings[i].push_back( ClusterInfo(featureInfo.at(cStart).getPosition(), clusterSize, molNumber) );
					} else {
						  // pick element at random
						size_t pos1 = rand() % clusterSize;
						  // first element
						selectedRings[i].push_back( ClusterInfo(featureInfo.at(cStart+pos1).getPosition(), clusterSize, molNumber) );
					}
				}
				  // go to next cluster
				cStart = cNext;
				++cID;
			}


			  // output histogram data to log
			*log <<"\n";
			*log <<" cluster size histogram :\n\tclusterSize\t#cluster\n";
			for (HIST::const_iterator c=clusterSizeHistogram.begin();c!=clusterSizeHistogram.end(); ++c) {
				*log <<"\t" <<c->first <<"\t" <<c->second <<"\n";
			}
			*log <<std::endl;
			*log <<" molecules per cluster histogram :\n\tmolPerCluster\t#cluster\n";
			for (HIST::const_iterator c=clusterMolHistogram.begin();c!=clusterMolHistogram.end(); ++c) {
				*log <<"\t" <<c->first <<"\t" <<c->second <<"\n";
			}
			*log <<std::endl;

			  // output selected SMILES to out stream
			for (size_t i=0; i<randomClusterLists; ++i) {
				*out <<"# cluster sample list " <<(i+1) <<"\n";
				for (std::vector<ClusterInfo>::const_iterator it=selectedRings.at(i).begin(); it!=selectedRings.at(i).end(); ++it) {

					  // get molecule/SMILES index of representative
					const size_t molID = trainSMILES.at(it->posInVector);

					  // get ring annotation within the representative
					annotateRing(molecules[molID], trainRing.at(it->posInVector), true);
					  // get ring annotated SMILES for output
					const std::string ringSMILES = SMILESwriter::getSMILES(molecules[molID]);
					  // undo ring annotation
					annotateRing(molecules[molID], trainRing.at(it->posInVector), false);

					  // output original SMILES
					*out <<SMILES.at(molID)
					  // output ring annotated SMILES
						<<" "
						<<ringSMILES
					  // output cluster size
						<<" "
						<<it->clusterSize
					  // output molecule number
						<<" "
						<<it->molNumber
					  // line break
						<<"\n";
				}
				*out <<std::endl;
			}


		}; break;
		case MODE_TEST :
		{

			/////////////////////////////////////////////////////////////
			// test mode
			/////////////////////////////////////////////////////////////

			  // check input stream
			if (boost::iequals(opts.getStrVal("smiles"),opts.getStrVal("model"))) {
				throw ArgException("cannot read both SMILES and model from same stream");
			}

			  // read model
			AP_NSPDK_Model model;
			model.modelID = modelID;

			  // check if default model to be used
			if (opts.getStrVal("model").size() == 1
				&& std::string("OoMmCc").find(opts.getStrVal("model").at(0)) != std::string::npos )
			{
				// load default model
				switch ( opts.getStrVal("model").at(0) ) {
				case 'O' :
				case 'o' :
					model = *(AP_NSPDK_Model::getInstance("OpenBabel:2013"));
					break;
				case 'M' :
				case 'm' :
					model = *(AP_NSPDK_Model::getInstance("Marvin:general:2013"));
					break;
				default :
					throw ArgException("unhandled default model character '"+opts.getStrVal("model")+"'");
					break;
				}
			} else {
				  // get model from stream

				  // set model input stream
				if (boost::iequals(opts.getStrVal("model"),"STDIN")) {
					  // read from STDIN
					inModel = &std::cin;
				} else {
					  // check if any model input given
					if (opts.getStrVal("model").size() == 0) {
						throw ArgException("test mode but no model input specified");
					}
					  // read from file
					inModelFile = new std::ifstream( opts.getStrVal("model").c_str()
												, std::ifstream::in );
					if (!inModelFile->is_open()) {
						std::ostringstream oss;
						oss	<<"cannot open model input file '" <<opts.getStrVal("model") <<"'";
						throw ArgException(oss.str());
					}
					inModel = inModelFile;
				}
				  // read model from stream
				readModel( *inModel, model );
			}
			if (opts.argExist("ringCentered") && !model.ringCentered) {
				throw ArgException("Parameter 'ringCentered' is not applicable for the used model");
			}

			  // create AP_NSPDK instance for feature generation using the model
			AP_NSPDK apNSPDK( & model );

			  // timing setup
			biu::Timer timer;
			double timeGenerating = 0.0;
			double timeTesting = 0.0;

			  // result counting
			enum Result { OVERALL, CORRECT
						, OVERALL_AROM, CORRECT_AROM
						, OVERALL_NON, CORRECT_NON
						, OVERALL_MOL, CORRECT_MOL
						, RES_NUM };
			std::vector< size_t > count( RES_NUM, 0 );

			  // test each ring of each molecule
			for (size_t i=0; i<SMILES.size(); ++i ) {

				count[OVERALL_MOL]++;
				  // create features of all rings within current molecule
				sgd::xvec_t testFeatures;
				sgd::yvec_t testTarget;
				timer.start();
				getFeatures( apNSPDK, SMILES.at(i), testFeatures, testTarget );
				timeGenerating += timer.stop();

				  // test each ring independently
				bool allRingsCorrect = true;
				for (size_t t=0; t<testFeatures.size(); ++t) {
					  // test instance
					timer.start();
					const double predictionValue = model.predictValue( testFeatures.at(t) );
					timeTesting += timer.stop();
					  // count instance
					count[OVERALL]++;
					  // report target and prediction if needed
					if (prediction != NULL) {
						*prediction <<testTarget.at(t) <<" " <<predictionValue <<"\n";
					}
					const bool isAromatic = predictionValue >= 0;
					if (testTarget.at(t) > 0) {
						count[OVERALL_AROM]++;
						if (isAromatic) {
							count[CORRECT]++;
							count[CORRECT_AROM]++;
						} else {
							allRingsCorrect = false;
						}
					} else {
						count[OVERALL_NON]++;
						if (!isAromatic) {
							count[CORRECT]++;
							count[CORRECT_NON]++;
						} else {
							allRingsCorrect = false;
						}
					}
				}

				  // count overall molecule prediction
				if (allRingsCorrect) {
					count[CORRECT_MOL]++;
				}
			}
			  // print statistics
			*out
				<<" time generating test   = " <<(size_t)timeGenerating <<" ms"
				<<"\n time testing           = " <<(size_t)timeTesting <<" ms"
				<<"\n overall rings correct  = " << count[CORRECT] << " / " <<count[OVERALL]
				<<" = " <<(double(count[CORRECT])/double(count[OVERALL])*100.0)
				<<" %"
				<<"\n aromatic rings correct = " << count[CORRECT_AROM] << " / " <<count[OVERALL_AROM]
				<<" = " <<(double(count[CORRECT_AROM])/double(count[OVERALL_AROM])*100.0)
				<<" %"
				<<"\n non-aromatic r correct = " << count[CORRECT_NON] << " / " <<count[OVERALL_NON]
				<<" = " <<(double(count[CORRECT_NON])/double(count[OVERALL_NON])*100.0)
				<<" %"
				<<"\n whole molecule correct = " << count[CORRECT_MOL] << " / " <<count[OVERALL_MOL]
				<<" = " <<(double(count[CORRECT_MOL])/double(count[OVERALL_MOL])*100.0)
				<<" %"
				<<"\n" <<std::endl;


		}; break;
		case MODE_APPLY :
		{

			/////////////////////////////////////////////////////////////
			// application mode
			/////////////////////////////////////////////////////////////

			  // check input stream
			if (boost::iequals(opts.getStrVal("smiles"),opts.getStrVal("model"))) {
				throw ArgException("cannot read both SMILES and model from same stream");
			}

			  // read model
			AP_NSPDK_Model model;
			model.modelID = modelID;

			  // check if default model to be used
			if (opts.getStrVal("model").size() == 1
				&& std::string("OoMmCc").find(opts.getStrVal("model").at(0)) != std::string::npos )
			{
				// load default model
				switch ( opts.getStrVal("model").at(0) ) {
				case 'O' :
				case 'o' :
					model = *(AP_NSPDK_Model::getInstance("OpenBabel:2013"));
					break;
				case 'M' :
				case 'm' :
					model = *(AP_NSPDK_Model::getInstance("Marvin:general:2013"));
					break;
				default :
					throw ArgException("unhandled default model character '"+opts.getStrVal("model")+"'");
					break;
				}
			} else {
				  // get model from stream

				  // set model input stream
				if (boost::iequals(opts.getStrVal("model"),"STDIN")) {
					  // read from STDIN
					inModel = &std::cin;
				} else {
					  // check if any model input given
					if (opts.getStrVal("model").size() == 0) {
						throw ArgException("test mode but no model input specified");
					}
					  // read from file
					inModelFile = new std::ifstream( opts.getStrVal("model").c_str()
												, std::ifstream::in );
					if (!inModelFile->is_open()) {
						std::ostringstream oss;
						oss	<<"cannot open model input file '" <<opts.getStrVal("model") <<"'";
						throw ArgException(oss.str());
					}
					inModel = inModelFile;
				}
				  // read model from stream
				readModel( *inModel, model );
			}

			if (opts.argExist("ringCentered") && !model.ringCentered) {
				throw ArgException("Parameter 'ringCentered' is not applicable for the used model");
			}

			  // create AP_NSPDK instance for feature generation using the model
			AP_NSPDK apNSPDK( & model );

			  // whether or not weighted GML output or SMILES is to be generated
			bool outputWeightGML = opts.argExist("outNodeWeight");
			bool multilineGML = false;

			  // run aromaticity perception for each molecule
			for (size_t i=0; i<SMILES.size(); ++i ) {

				  // get molecule graph
				std::pair< Molecule, int > parse = SMILESparser::parseSMILES( SMILES.at(i) );

				  // proton handling
				Molecule mol;
				MoleculeUtil::copy( parse.first, mol );
				MoleculeUtil::fillProtons( mol );

				  // apply aromaticity perception
				apNSPDK.correctAromaticity( mol, false );

				if (outputWeightGML) {
					  // write weighted GML output
					printWeightedGML( *out, apNSPDK, mol, multilineGML );
				} else {
					  // write SMILES output
					*out <<ggl::chem::SMILESwriter::getSMILES( mol );
				}
				  // line break after each molecule
				*out <<std::endl;
			}

		}; break;
		case MODE_UNDEF : {
			throw ArgException("no mode selected"); // should never occur
		}; break;
		default:
			throw ArgException("unknown mode selected"); // should never occur
		}

	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// final stream handling
	//////////////////////////////////////////////////////////////

	inSMILES = &std::cin;
	inModel = &std::cin;
	out = &std::cout;
	log = &std::cerr;
	prediction = NULL;
	if (inSMILESFile != NULL)	{ inSMILESFile->close(); delete inSMILESFile; }
	if (inModelFile != NULL)	{ inModelFile->close(); delete inModelFile; }
	if (outFile != NULL)		{ outFile->close(); delete outFile; }
	if (logFile != NULL)		{ logFile->close(); delete logFile; }
	if (predictionFile != NULL)	{ predictionFile->close(); delete predictionFile; }

	return exitValue;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

/*!
 * Initialises allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	infoText = "\n"
		"Reads a list of SMILES and trains/tests/applies an aromaticity"
		" model using NSPDK and the SGD framework.\n"
		"\n"
		;

	allowedArgs.push_back(biu::COption(
							"mode", false, biu::COption::STRING,
							"Running mode : 'train' : train a model using the given parameters;"
							" 'xval' : does an 'xval'-cross validation using the given model training parameters and prints target-prediction value pairs to 'out';"
							" 'cluster' : clusters the training data based on feature identity for a model using the given parameters;"
							" 'test' : tests the prediction of ring wise and for whole molecules (model needed);"
							" 'apply' : applies the given model to each molecule and prints the corrected SMILES (model needed)"
							));
	allowedArgs.push_back(biu::COption(
							"smiles", true, biu::COption::STRING,
							"SMILES input file name or 'STDIN' when to read from standard input",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"out", true, biu::COption::STRING,
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"log", true, biu::COption::STRING,
							"Log file name or 'STDERR'/'STDOUT' when to write to standard error / output",
							"STDERR"));
	allowedArgs.push_back(biu::COption(
							"prediction", true, biu::COption::STRING,
							"Prediction file name or 'STDERR'/'STDOUT' when to write to standard error / output"));
	allowedArgs.push_back(biu::COption(
							"model", true, biu::COption::STRING,
							"models input file name, or 'STDIN' when to read from standard input, or a single letter code for one of the default models : (O)penBabel or (M)arvin general mode.",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"modelID", true, biu::COption::STRING,
							"the models ID to set",
							"NSPDK-SGD-aromaticity-model"));
	allowedArgs.push_back(biu::COption(
							"maxRingSize", true, biu::COption::INT,
							"the maximal ring size to be considered for aromatic ring training/testing",
							"13"));
	allowedArgs.push_back(biu::COption(
							"ringCentered", true, biu::COption::BOOL,
							"NSPDK: if given, the NSPDK features have to touch the ring of interest in at least one node"));
	allowedArgs.push_back(biu::COption(
							"nspdkMaxRingDistance", true, biu::COption::INT,
							"NSPDK: IF(ringCentered) : the distance from ring nodes to be considered for feature generation"));
	allowedArgs.push_back(biu::COption(
							"nspdkMaxDistance", true, biu::COption::INT,
							"NSPDK: the maximal distance between nodes to be considered for feature generation",
							"4"));
	allowedArgs.push_back(biu::COption(
							"nspdkMaxRadius", true, biu::COption::INT,
							"NSPDK: the maximal radius around nodes to be considered for feature generation",
							"2"));
	allowedArgs.push_back(biu::COption(
							"nspdkFeatureBitSize", true, biu::COption::INT,
							"NSPDK: the number of bits used to encode the NSPDK features",
							"22"));
	allowedArgs.push_back(biu::COption(
							"lambda", true, biu::COption::DOUBLE,
							"[mode=train] the lambda factor used by the SGD framework",
							"1e-4"));
	allowedArgs.push_back(biu::COption(
							"epochs", true, biu::COption::INT,
							"[mode=train] the number of training epochs performed by the SGD framework",
							"5"));
	allowedArgs.push_back(biu::COption(
							"outClustLists", true, biu::COption::INT,
							"[mode=cluster] number of random cluster sample lists to be reported",
							"1"));
	allowedArgs.push_back(biu::COption(
							"xval", true, biu::COption::INT,
							"[mode=xval] sets the X for X-cross validation, ie. the number of groups to randomly split the input into",
							"10"));
	allowedArgs.push_back(biu::COption(
							"outNodeWeight", true, biu::COption::BOOL,
							"[mode=apply] output is given in GML format where for each ring the normalized feature weights for view point nodes are provided"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////

void getFeatures( AP_NSPDK & ap
				, const std::string & molSMILES
				, sgd::xvec_t& trainFeatures
				, sgd::yvec_t& trainTarget )
{

	  // get molecule graph
	std::pair< Molecule, int > parse = SMILESparser::parseSMILES( molSMILES );

	  // fill protons
	Molecule mol;
	MoleculeUtil::copy( parse.first, mol );
	MoleculeUtil::fillProtons( mol );

	  // call final function
	getFeatures(ap, mol, trainFeatures, trainTarget, NULL);

}


//////////////////////////////////////////////////////////////////////////

void getFeatures( AP_NSPDK & ap
				, const Molecule & mol
				, sgd::xvec_t& trainFeatures
				, sgd::yvec_t& trainTarget
				, std::vector< sgm::RingReporter::RingList > * trainRing
				)
{

	  // get features
	std::vector< std::pair< sgm::RingReporter::RingList, nspdk::SVector> > features = ap.getFeatures( mol );

	  // create a dummy graph to get access to the ring labels
	ggl::chem::Molecule_Graph molGraph(mol);

	  // feed all rings to the feature vector
	for (size_t i=0; i<features.size(); ++i ) {

		  // check if aromatic ring via bond labels
		sgm::RingReporter::RingList::const_iterator cur = features.at(i).first.begin();
		sgm::RingReporter::RingList::const_iterator last = features.at(i).first.begin();
		bool isAromatic = true;
		for (++cur; isAromatic && cur!=features.at(i).first.end();++cur,++last) {
			assert( molGraph.getEdgesBegin(*last,*cur) != molGraph.getEdgesEnd(*last,*cur) /*no edge present*/);
			isAromatic = MoleculeUtil::getBondData(molGraph.getEdgesBegin(*last,*cur)->getEdgeLabel())->isAromatic != 0;
		}

		  // store features
		trainFeatures.push_back( features.at(i).second );
		  // store training target
		trainTarget.push_back( isAromatic ? +1 : -1 );
		  // store ring information
		if (trainRing != NULL) {
			trainRing->push_back( features.at(i).first );
		}


	}
}

//////////////////////////////////////////////////////////////////////////

/*!
 * Parses SMILES strings from stream. Each line should contain only ONE SMILES
 * string. Leading and tailing whitespaces are ignored.
 *
 * @param in the stream to read from
 * @param toFill the STL inserter to add the found SMILES to
 *
 */
void
parseSMILES(	std::istream & in
				, SMILES_inserter & toFill
				, const size_t linesRead ) throw(std::exception)
{
	std::string line;
	size_t lineNumber = linesRead;
	const std::string WHITESPACES = std::string(" \t\n");
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;

		  // check if line is not empty
		if (line.size() > 0) {
			size_t start = line.find_first_not_of(WHITESPACES);
			  // check if line is not only filled with blanks
			if (start != std::string::npos ) {
				size_t end = line.find_last_not_of(WHITESPACES);
				  // add SMILES string to inserter
				*toFill = line.substr( start, end - start + 1 );
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////

void
readModel( std::istream & in, AP_NSPDK_Model & model )
{
	in
		>> model.ringCentered
		>> model.nspdk_maxRingDistance
		>> model.nspdk_maxDistance
		>> model.nspdk_maxRadius
		>> model.nspdk_featureBitSize
		>> model.maxRingSize
		>> model.bias
		>> model.wscale
		>> model.w
	;

	model.ringEdgeLabel = "*~";
	model.ringNodeLabelPrefix = "*";
	model.ringViewLabelPrefix = "?";
}

//////////////////////////////////////////////////////////////////////////

void
printModel( std::ostream & out, const AP_NSPDK_Model & model )
{
	out
		<< model.ringCentered <<" "
		<< model.nspdk_maxRingDistance <<" "
		<< model.nspdk_maxDistance <<" "
		<< model.nspdk_maxRadius <<" "
		<< model.nspdk_featureBitSize <<" "
		<< model.maxRingSize <<" "
		<< model.bias <<" "
		<< model.wscale <<" "
		<< model.w <<"\n"
	;
}

//////////////////////////////////////////////////////////////////////////

void
annotateRing( Molecule & mol, const sgm::RingReporter::RingList & ring, bool addAnnotation )
{

	const std::string ringClass = ":1";

	  // label access
	boost::property_map< ggl::chem::Molecule, ggl::chem::PropNodeLabel >::type
		nodeLabel = boost::get( ggl::chem::PropNodeLabel(), mol );

	  // update ring vertices
	sgm::RingReporter::RingList::const_iterator cur = ring.begin();
	for (++cur; cur!=ring.end();++cur) {

		  // get vertex
		Molecule::vertex_descriptor node = boost::vertex( *cur, mol );

		  // get current atom label
		std::string atomLabel = nodeLabel[node];

		  // alter the atom label accordingly
		if (addAnnotation) {
			  // ensure the ring class is not present
			assert( atomLabel.find(ringClass) == std::string::npos );
			atomLabel += ringClass;
		} else {
			  // ensure the ring class is present
			assert( atomLabel.find(ringClass) != std::string::npos );
			atomLabel = atomLabel.substr(0,atomLabel.find(ringClass));

		}

		  // set altered label
		nodeLabel[ node ] = atomLabel;
	}

}


//////////////////////////////////////////////////////////////////////////

void
printWeightedGML( std::ostream & out, AP_NSPDK & apNSPDK, const Molecule & mol, const bool multilineGML )
{

	std::string linesep = (multilineGML ? "\n" : "");
	std::string indent = (multilineGML ? "  " : "");

	  // get weights for each ring
	std::vector< std::pair< AP_NSPDK::RingDescriptor, AP_NSPDK::RingWeightVec> >
		weightsPerRing = apNSPDK.getNodeWeights( mol );

	  // generate additional GML information for each ring
	std::stringstream ringGmlEncoding;
	  // change bond and node labels
	  // for each predicted ring
	for (size_t i=0; i<weightsPerRing.size(); ++i)
	{
		const AP_NSPDK::RingDescriptor & ring = weightsPerRing.at(i).first;
		const AP_NSPDK::RingWeightVec & ringWeights = weightsPerRing.at(i).second;

		  // general ring informationring
		ringGmlEncoding <<"ring [" <<linesep
						<<indent <<"label \""
						<<(ring.predState == AP_NSPDK::RingDescriptor::Aromatic ? "aromatic" : "non-aromatic" )
						<<":"
						<<ring.predCertainty
						<<"\""<<linesep
		  // ring nodes
						<<" ringNodes [ "
						;
		  // start from second node, since first and last node are the same
		sgm::RingReporter::RingList::const_iterator ringNode = ring.ring.begin();
		for (ringNode++;ringNode != ring.ring.end(); ++ringNode )
		{
			ringGmlEncoding <<"id " <<*ringNode <<" ";
		}
		ringGmlEncoding	<<"]" <<linesep;

		  // print node weights
		for (AP_NSPDK::RingWeightVec::const_iterator weight = ringWeights.begin();
				weight != ringWeights.end(); ++weight )
		{
			ringGmlEncoding <<indent <<"nodeWeight [ id " <<weight->first <<" label \"" <<weight->second <<"\" ]" <<linesep;
		}

		  // print GML end
		ringGmlEncoding <<"]" <<linesep;

	}

	  // write final molecule to stream in GML format
	std::string additionalGML = ringGmlEncoding.str();
	Graph_GML_writer::write( out, mol, multilineGML, &additionalGML );

}

//////////////////////////////////////////////////////////////////////////

bool
SVector_less(const sgd::SVector& v, const sgd::SVector & vCompare)
{
		  // get iterators over sparse entries
		const sgd::SVector::Pair *p1 = v, *p2 = vCompare;

		// skip all equal elements
		for( ; p1->i >= 0 && p2->i >= 0 && p1->i == p2->i && p1->v == p2->v; p1++,p2++)
		{ /* just go to next entries */}

		  // check first non-equal element
		if (p1->i < p2->i ) {
			  // if (p1->i < 0) one additional larger entry in vCompare,
			  // else one entry in v existing that is 0 in vCompare
			return (p1->i < 0 && p2->v > 0) || (p1->i > 0 && p1->v < 0);
		} else if (p1->i > p2->i) {
			  // if (p2->i < 0) one additional larger entry in v,
			  // else one entry in vCompare existing that is 0 in v
			return (p2->i < 0 && p1->v < 0) || (p2->i > 0 && p2->v > 0);
		} else {
			// p1->i == p2->i
			if (p1->i < 0) {
				  // all the same
				return false;
			} else {
				  // check if v is smaller
				return p1->v < p2->v;
			}
		}

}


//////////////////////////////////////////////////////////////////////////
