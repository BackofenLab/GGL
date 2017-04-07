

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <exception>
#include <iterator>

#include <boost/algorithm/string.hpp>

#include "biu/OptionParser.hh"

#include <sgm/GM_vf2.hh>
#include <sgm/MR_Storing.hh>
#include <sgm/SubGraph.hh>
#include "sgm/RP_Hanser96.hh"

#include <ggl/Graph.hh>
#include <ggl/Graph_GML_writer.hh>
#include <ggl/chem/Molecule.hh>
#include <ggl/chem/SMILESwriter.hh>
#include <ggl/chem/SMILESparser.hh>


#include "version.hh"
#include "toyChemUtil.hh"

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


void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );

//////////////////////////////////////////////////////////////////////////

/*!
 *
 * Dummy aromaticity perception to access rings to predict
 *
 */
class AP_enumerate
	: public AromaticityPerception
{

protected:

	  //! defines an edge; NOTE: use as ordered pairs, i.e. first <= second
	typedef AromaticityPerception::Edge AromaticEdge;
	  //! defines an aromatic ring via the set of aromatic edges
	typedef AromaticityPerception::EdgeSet AromaticEdgeSet;

public:
	  //! ring describing class
	typedef AromaticityPerception::RingDescriptor RingDescriptor;

protected:
	  //! the aromatic ring container to fill
	using AromaticityPerception::aromaticEdges;


	  //! container that will hold all rings of the current molecule
	  //! -> this list is later pruned to the rings of interest
	using AromaticityPerception::allRings;

	const size_t maxRingSize;

public:

	AP_enumerate(const size_t maxRingSize)
		: maxRingSize(maxRingSize)
	{}

	AP_enumerate( const AP_enumerate & toCopy )
		: AromaticityPerception(toCopy)
		, maxRingSize(toCopy.maxRingSize)
	{}

	virtual ~AP_enumerate()
	{}

	const std::vector< RingDescriptor* > &
	enumerateRings( const Molecule & mol ) {
		  // clear temporary data
		clearData();
		  // identify all rings
		findAllRings( mol );

		  // prune rings to remove fused representation
		pruneFusedRings( allRings );

		  // prune rings that cannot be aromatic
		pruneNonSingleDoubleBondRings( allRings, mol );

		return allRings;
	}

	virtual
	void
	findAllRings( const Molecule & mol ) {

		  // store all rings up to the maximally predictable size
		sgm::RP_Hanser96 ringFinder;
		  // setup graph interface
		Molecule_Graph molGraph(mol);
		  // run ring perception
		ringFinder.findRings( molGraph, *this, maxRingSize );

	}

	void
	identifyAromaticEdges( const Molecule & mol )
	{}


	virtual
	AP_enumerate *
	clone() const
	{
		return new AP_enumerate(*this);
	}

protected:  // member functions

	  /*!
	   * Computes the NSPDK feature vector for a given ring within a
	   * molecule. Note, the molecule graph has to show already the
	   * according ring relabeling.
	   * @param molGraph the molecule graph containing the ring
	   * @param ring the rings to be described
	   * @param aromaticityModel the model to use
	   */
	static
	nspdk::SVector
	getFeaturesOfRing( nspdk::GraphClass &molGraph
						, const RingDescriptor & ring
						, const AP_NSPDK_Model& aromaticityModel);

	  /*!
	   * Does an aromaticity prediction for the given ring.
	   * @param ring the ring to predict
	   * @return true if the ring is aromatic; false otherwise
	   */
	bool
	isAromaticRing( const RingDescriptor & ring );


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
};



//////////////////////////////////////////////////////////////////////////


int main( int argc, char** argv ) {

	//////////////////////////////////////////////////////////////
	// variables
	//////////////////////////////////////////////////////////////

	std::istream* in = NULL;
	std::ifstream* inFile = NULL;
	size_t inLine = 0;

	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;

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
		if (boost::iequals(opts.getStrVal("graph"),"STDIN")) {
			  // read from STDIN
			in = &std::cin;
		} else {
			  // read from file
			inFile = new std::ifstream( opts.getStrVal("graph").c_str()
										, std::ifstream::in );
			if (!inFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open input file '" <<opts.getStrVal("graph") <<"'";
				throw ArgException(oss.str());
			}
			in = inFile;
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

		const int maxRingSize = opts.getIntVal("maxRingSize");
		if (maxRingSize < 3) {
			throw ArgException("'maxRingSize' has to be >= 3");
		}

		//////////////////////////////////////////////////////////////
		// parse molecules from input
		//////////////////////////////////////////////////////////////

		  // get master molecule to match onto
		SMILES_container parsedMols;
		ggl::chem::GroupMap dummyGroups;
		parseMolGML( *in, parsedMols, inLine, dummyGroups );
		if (parsedMols.size() != 1)  {
			throw ArgException("could not parse only one molecule from 'graph'");
		}
		Molecule mol = *(parsedMols.begin()->second);


		//////////////////////   RING ENUMERATION  /////////////////////////////////////

		AP_enumerate ringEnumerator(maxRingSize);
		const std::vector< AP_enumerate::RingDescriptor* > & allRings = ringEnumerator.enumerateRings( mol );


		//////////////////////   OUTPUT  /////////////////////////////////////

		for (size_t i=0; i< allRings.size(); ++i) {

			const sgm::RingReporter::RingList & nodes = allRings.at(i)->ring;

			for (sgm::RingReporter::RingList::const_iterator n = nodes.begin(); n!= nodes.end(); ++n) {
				*out <<" "<<*n;
			}
			*out <<std::endl;

		}


	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// final stream handling
	//////////////////////////////////////////////////////////////

	in = &std::cin;
	out = &std::cout;
	if (inFile != NULL)		{ inFile->close(); delete inFile; }
	if (outFile != NULL)	{ outFile->close(); delete outFile; }

	return exitValue;
}



//////////////////////////////////////////////////////////////////////////

/*!
 * Initialises allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	infoText = "\n"
		"Enumerates all aromatic ring candidates for a given molecule.\n"
		;

	allowedArgs.push_back(biu::COption(
							"graph", true, biu::COption::STRING,
							"The name of the file that holds the molecule encoding in GML or 'STDIN' when to read GML from standard input",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"maxRingSize", true, biu::COption::INT,
							"The maximal ring size to be considered",
							"15"));
	allowedArgs.push_back(biu::COption(
							"out", true, biu::COption::STRING,
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
