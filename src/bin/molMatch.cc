

#include <iostream>
#include <fstream>
#include <sstream>
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


class NonAromNodeGraph : public Molecule_Graph {
protected:

	  //! non-aromatic node labels for faster access
	std::vector< std::string > nodeLabels;

public:

	  //! Construction of the interface graph.
	  //! @param graph the boost graph to use as internal data structure
	NonAromNodeGraph( const Molecule_Graph::InternalBoostGraph & graph )
		: Molecule_Graph(graph)
		, nodeLabels()
	{
		  // fill nodeLabels
		nodeLabels.resize(getNodeNumber());
		for (size_t i=0; i<getNodeNumber(); ++i) {
			  // get super class result
			std::string curLabel = Molecule_Graph::getNodeLabel(i);

			const MoleculeUtil::AtomLabelData * atomData = MoleculeUtil::getAtomData( curLabel );
			  // check if aromatic
			if (atomData->isAromatic == 0) {
				  // non-aromatic -> return unchanged
				nodeLabels[i] = curLabel;
			} else {
				  // aromatic -> create new label
				std::string atomLabel = MoleculeUtil::getAtom( curLabel );
				assert( MoleculeUtil::getAromaticPendant( atomLabel )!= NULL );
				atomLabel = *(MoleculeUtil::getAromaticPendant( atomLabel ));
				  // create new according atom label and store
				nodeLabels[i] =  MoleculeUtil::getComplexAtomLabel( atomLabel
										, MoleculeUtil::getProtons(curLabel)
										, MoleculeUtil::getCharge(curLabel)
										, MoleculeUtil::getClass(curLabel) );
			}
		}
	}


	  //! Access to the non-aromatic label of a specified node
	  //! @param i the index of the node of interest
	  //! @return a string representation of the node label
	virtual
	std::string
	getNodeLabel(const IndexType & i) const {
		assert( i < nodeLabels.size() );
		return nodeLabels.at(i);
	}

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

		const int maxMappings = opts.getIntVal("max");
		if (maxMappings < 1) {
			throw ArgException("'max' has to be >= 1");
		}

		// check outMode parameter setup : non-empty and within allowed chars
		const std::string allowedModeChars = "gGlL";
		const char outMode = opts.getCharVal("outMode");
		if (allowedModeChars.find(outMode) == std::string::npos) {
			std::ostringstream oss;
			oss	<<"outMode '" <<opts.getCharVal("outMode") <<"' is not supported";
			throw ArgException(oss.str());
		}
		  // check for single line GML
		const bool oneLineGML = std::string("lL").find(outMode) != std::string::npos;

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
		Molecule molMaster = *(parsedMols.begin()->second);

		  // get SMILES of molecule to match
		std::string mol2matchSMILES = opts.getStrVal("smiles");
		if (mol2matchSMILES.empty()) {
			throw ArgException("given SMILES is empty");
		}
		  // convert to molecule
		std::pair<Molecule,int> ret = ggl::chem::SMILESparser::parseSMILES( mol2matchSMILES );
		if (ret.second >= 0) {
			std::stringstream error;
			error <<"cannot parse SMILES at position " <<ret.second <<":\n"
					<<mol2matchSMILES.substr(0,ret.second) <<"\n";
			for (int i=0; i<ret.second; ++i) { error <<" "; }
			error <<mol2matchSMILES.substr(ret.second) <<"\n";
			throw ArgException(error.str());
		}
		Molecule mol2match = ret.first;

		  // fill protons if requested
		if (!opts.argExist("noProtonFilling")) {
			MoleculeUtil::fillProtons( mol2match );
		}

		//////////////////////   MATCHING  /////////////////////////////////////

		  // get graph interface representations that ignores aromaticity on node labels
		NonAromNodeGraph graphMaster( molMaster ), graph2match( mol2match );

		if (opts.argExist("v")) {
			std::cerr <<"\n master graph : \n" <<graphMaster <<"\n pattern graph : \n" <<graph2match <<"\n"; // DEBUG OUTPUT
		}

		  // setup storage for all mappings
		sgm::MR_Storing::Storage storage;
		sgm::MR_Storing matchStoring(storage);

		  // get pattern that ignores aromatic edges during matching
		std::string aromaticEdge = ":";
		sgm::Pattern graphPattern = sgm::Pattern(graph2match, aromaticEdge);

		  // compute mapping
		sgm::GM_vf2 matcher;
		matcher.findMatches( graphPattern, graphMaster, matchStoring, maxMappings);

		//////////////////////   OUTPUT  /////////////////////////////////////

		Molecule_Graph molGraph2match(mol2match);
		for (size_t i=0; i< storage.size(); ++i) {

			sgm::Match atomMapping = storage.at(i);
			sgm::SubGraph mappedGraph( molGraph2match, atomMapping );

			Molecule mappedMol;
			MoleculeUtil::copy( mappedGraph, mappedMol );

			  // write mapped molecule in GML to stream
			ggl::Graph_GML_writer::write( *out, mappedMol, !oneLineGML );
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
		"Computes a mapping of two molecule representations ignoring aromaticity.\n"
		;

	allowedArgs.push_back(biu::COption(
							"graph", false, biu::COption::STRING,
							"The name of the file that holds the 'master' molecule encoding in GML or 'STDIN' when to read GML from standard input",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"smiles", false, biu::COption::STRING,
							"The encoding molecule to map in SMILES format"
							));
	allowedArgs.push_back(biu::COption(
							"out", true, biu::COption::STRING,
							"Output file name where to write the GML encoding or 'STDOUT' when to write GML to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"outMode", true, biu::COption::CHAR,
							"Output format : (G)ML of one (L)ine GML",
							"G"));
	allowedArgs.push_back(biu::COption(
							"max", true, biu::COption::INT,
							"The maximal number of mappings to output",
							"1"));
	allowedArgs.push_back(biu::COption(
							"noProtonFilling", true, biu::COption::BOOL,
							"Dont add protons to the molecule given in SMILES before matching it"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(
							"v", true, biu::COption::BOOL,
							"Enables verbose output to STDERR"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
