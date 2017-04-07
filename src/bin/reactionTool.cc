/*
 * Enables the conversion of a GML molecule graph description into a SMILES
 * string representation and vice versa.
 * Furthermore, enables the correction and check of the molecule.
 *
 *  Created on: 22.03.2012
 *      Author: mmann
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <exception>
#include <iterator>

#include <boost/algorithm/string.hpp>

#include "biu/OptionParser.hh"

#include <sgm/SubGraph.hh>

#include <ggl/Rule_GML_writer.hh>
#include <ggl/RuleGraph.hh>

#include <ggl/chem/ChemRule.hh>
#include <ggl/chem/SMILESparser.hh>
#include <ggl/chem/ClassIDGraph.hh>

#include "version.hh"
#include "toyChemUtil.hh"

	using namespace ggl;
	using namespace ggl::chem;


typedef std::vector< ggl::chem::ChemRule > RULE_container;


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

/*
 * Parses SMILES strings from stream. Each line should contain only ONE SMILES
 * string. Leading and tailing whitespaces are ignored.
 *
 * @param in the stream to read from
 * @param toFill the inserter to add the found SMILES to
 *
 */
void
parseReactionSMILES(	std::istream & in
						, RULE_container & toFill
						, const size_t linesRead
						, const bool pruneClassID = true ) throw(std::exception);

void
parseReactionSMILES(	const std::string & inSource
					, std::vector<ggl::chem::ChemRule> & toFill
					, const bool pruneClassID = true ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

/**
 * writes all components of a given multi-molecule graph as SMILES to stream
 * while molecule are separated by a '.'.
 *
 * @param multiMol the multi-molecule graph to print
 * @param out the output stream to write to
 * @param ignoreProtons whether or not protons are to be ignored in SMILES
 *        output or not
 * @param groups the groups that will be compressed if found
 */
void
writeMoleculesToStream(	const sgm::Graph_Interface & multiMol
						, std::ostream & out
						, const bool ignoreProtons
						, const ggl::chem::GroupMap & groups );

//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////


int main( int argc, char** argv ) {

	//////////////////////////////////////////////////////////////
	// variables
	//////////////////////////////////////////////////////////////

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


		  // whether or not protons are to be ignored in SMILES output
		const bool ignoreProtons = !opts.argExist("noProtonPruning");

		const bool noClassPruning = opts.argExist("noClassPruning");


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


		// check outMode parameter setup : non-empty and within allowed chars
		const std::string allowedModeChars = "gGsScC";
		const char outMode = opts.getCharVal("outMode");
		if (allowedModeChars.find(outMode) == std::string::npos) {
			std::ostringstream oss;
			oss	<<"outMode '" <<opts.getCharVal("outMode") <<"' is not supported";
			throw ArgException(oss.str());
		}


		//////////////////////////////////////////////////////////////
		// parse groups from input
		//////////////////////////////////////////////////////////////

		ggl::chem::GroupMap moleculeGroups;

		if (opts.argExist("groups")) {
			if (opts.getStrVal("groups") == opts.getStrVal("in")) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and input from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}
			  // parse molecule group definitions that can be abbreviated in
			  // rules and molecules
			parseGroups( opts.getStrVal("groups"), moleculeGroups );
//			printGroups( std::cerr, moleculeGroups );
		}
		  // check wether or not groups are to be inserted into the output
		ggl::chem::GroupMap noGroups;
		ggl::chem::GroupMap * parseGroups = &moleculeGroups;
		ggl::chem::GroupMap * outGroups = &noGroups;
		if (opts.argExist("noGroupInsertion")) {
			parseGroups = &noGroups;
			outGroups = &moleculeGroups;
		}


		RULE_container parsedRules;
		if (opts.getCharVal("inMode") == 'S' || opts.getCharVal("inMode") == 's') {
			parseReactionSMILES( opts.getStrVal("in"), parsedRules, !noClassPruning );
		} else if (opts.getCharVal("inMode") == 'G' || opts.getCharVal("inMode") == 'g') {
			parseRules( opts.getStrVal("in"), parsedRules, *parseGroups );
		} else if (opts.getCharVal("inMode") == 'C' || opts.getCharVal("inMode") == 'c') {
			parseRules( opts.getStrVal("in"), parsedRules, *parseGroups, true );
		} else {
			std::ostringstream oss;
			oss	<<"Input format '" <<opts.getStrVal("inMode") <<"' is unknown";
			throw ArgException(oss.str());
		}

		  // do sanity checks
		if (!opts.getBoolVal("noInputCheck")) {
			for (size_t i=0; i<parsedRules.size() ; ++i ) {
				size_t conStatus = parsedRules.at(i).isConsistent();
				if( conStatus != ggl::chem::ChemRule::C_Consistent ) {
					(*out) <<"\n PROBLEM : reaction "
							<<i <<" is not chemically correct"
								" or contains unsupported properties:\n";
					parsedRules.at(i).decodeConsistencyStatus( conStatus, (*out) );
					(*out) <<std::endl;
				}
			}
		}



		//////////////////////   OUTPUT  /////////////////////////////////////

		  // write GML results to stream
		switch( opts.getCharVal("outMode") ) {

		case 'G' : case 'g' :
			// write in GML format
			for( RULE_container::const_iterator rule=parsedRules.begin();
					rule != parsedRules.end(); ++rule )
			{
				  // write molecule in GML to stream
				ggl::Rule_GML_writer::write( *out, *(rule), true );
				*out <<std::endl;
			}
			break;

		case 'S' : case 's' :
			// write SMIRKS results to stream
			for( RULE_container::const_iterator rule=parsedRules.begin();
					rule != parsedRules.end(); ++rule )
			{
				  // the SMIRKS representation to fill
				std::stringstream smirks;
				  // get the educt molecules
				ggl::chem::LeftSidePattern educts(*rule);
//				if (noClassPruning) {
					  // never prune class IDs in SMIRK output, otherwise reaction is not represented correctly
					ClassIDGraph classIDwrapper( educts, 2 );
					writeMoleculesToStream( classIDwrapper, smirks, ignoreProtons, *outGroups );
//				} else {
//					writeMoleculesToStream( educts, smirks, ignoreProtons, *outGroups );
//				}

				  // append separators
				smirks << ">>";

				  // get the product molecules
				ggl::chem::RightSidePattern products(*rule);
				if (noClassPruning) {
					ClassIDGraph classIDwrapper( products, 2 );
					writeMoleculesToStream( classIDwrapper, smirks, ignoreProtons, *outGroups );
				} else {
					writeMoleculesToStream( products, smirks, ignoreProtons, *outGroups );
				}

				  // write SMIRKS to stream
				*out <<smirks.str() <<" " <<rule->getID() <<std::endl;
			}
			break;

		case 'C' : case 'c' :
			// write in Compact GML format
			for( RULE_container::const_iterator rule=parsedRules.begin();
					rule != parsedRules.end(); ++rule )
			{
				  // write molecule in GML to stream
				ggl::Rule_GML_writer::writeCompact( *out, *(rule), true );
				*out <<std::endl;
			}
			break;

		default:
			// unknown mode
			std::ostringstream oss;
			oss	<<"outMode '" <<opts.getCharVal("outMode") <<"'"
				<<" is not supported";
			throw ArgException(oss.str());

		}

	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// final stream handling
	//////////////////////////////////////////////////////////////

	out = &std::cout;
	if (outFile != NULL)		{ outFile->close(); delete outFile; }

	return exitValue;
}



//////////////////////////////////////////////////////////////////////////

/*!
 * Initializes allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	infoText = "\n"
		"Enables the check and conversion of chemical reactions in different formats link SMILES or GML.\n"
		;

	allowedArgs.push_back(biu::COption(
							"in", true, biu::COption::STRING,
							"Input file name or 'STDIN' when to read from standard input",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"inMode", true, biu::COption::CHAR,
							"Input format : (S)MIRKS reaction SMILES, (G)ML, or (C)ompacted GML",
							"S"));
	allowedArgs.push_back(biu::COption(
							"groups", true, biu::COption::STRING,
							"Predefined molecule groups abbreviated in rules and molecules in GML format : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(
							"out", true, biu::COption::STRING,
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"outMode", true, biu::COption::CHAR,
							"Output format : (G)ML, (S)MIRK reaction SMILES, or (C)ompacted GML",
							"G"));
	allowedArgs.push_back(biu::COption(
							"noInputCheck", true, biu::COption::BOOL,
							"Dont check the input rule for consistency (atom/bond label, ...)"));
	allowedArgs.push_back(biu::COption(
							"noClassPruning", true, biu::COption::BOOL,
							"Dont prune the class ID information from the atom labels (useful for debug purpose)"));
	allowedArgs.push_back(biu::COption(
							"noProtonPruning", true, biu::COption::BOOL,
							"SMILES output : Dont prune unnecessary hydrogens"));
	allowedArgs.push_back(biu::COption(
							"noGroupInsertion", true, biu::COption::BOOL,
							"Dont replace molecule groups with their identifiers (from 'groups' parameter) before producing the output"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"Version information"));
}


//////////////////////////////////////////////////////////////////////////


/*
 * Parses SMILES strings from stream. Each line should contain only ONE SMILES
 * string. Leading and tailing whitespaces are ignored.
 *
 * @param in the stream to read from
 * @param toFill the inserter to add the found SMILES to
 *
 */
void
parseReactionSMILES(	std::istream & in
						, RULE_container & toFill
						, const size_t linesRead
						, const bool pruneClassID ) throw(std::exception)
{
	std::string line;
	std::string smiles;
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
				  // get SMILES string
				smiles = line.substr( start, end - start + 1 );
				  // check if WHITESPACE within SMILES
				if (smiles.find(WHITESPACES) != std::string::npos) {
					std::ostringstream oss;
					oss	<<"parsing error in line " <<(linesRead+1)
						<<" in SMILES string '"
						<<smiles
						<<"' contains whitespaces";
					throw ArgException(oss.str());
				}
				  // parse SMILES to graph
				ggl::chem::ChemRule::CoreGraph result;
				try {
					result = ggl::chem::SMILESparser::parseReactionSMILES(smiles,pruneClassID);
				} catch (std::exception& e) {
					  // error handling
					std::ostringstream oss;
					oss	<<"parsing error in line " <<(linesRead+1)
						<<" in SMILES string '"
						<<smiles
						<<"' : " <<e.what();
					throw ArgException(oss.str());
				}

				  // store parsed rule graph
				toFill.push_back(ChemRule(result,"SMILES rule"));
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////


void
parseReactionSMILES(	const std::string & inSource
					, std::vector<ggl::chem::ChemRule> & toFill
					, const bool pruneClassID ) throw(std::exception)
{
	if (inSource.compare("STDIN") == 0) {
		size_t linesRead = 0;
		parseReactionSMILES( std::cin, toFill, linesRead, pruneClassID );
		return;
	}

	  // split inSource into file name list if any
	std::vector< std::string > files = splitToFileNames( inSource );

	  // read all files
	size_t i=0;
	for (	std::vector< std::string >::const_iterator file=files.begin();
			file != files.end(); ++file)
	{
		i++;
		  // ensure there is no standard stream within the list
		if (file->compare("STDIN") == 0) {
			std::ostringstream oss;
			oss	<<"RULES input list : the " <<i
				<<". entry is 'STDIN', but it has to be a file list!";
			throw ArgException(oss.str());
		}
		 // open current file
		std::ifstream inFile( file->c_str(), std::ifstream::in );
		if (!inFile.is_open()) {
			std::ostringstream oss;
			oss	<<"cannot open " <<i <<". RULES input file '" <<(*file) <<"'";
			throw ArgException(oss.str());
		}
		 // parse current file
		size_t linesRead = 0;
		parseReactionSMILES( inFile, toFill, linesRead, pruneClassID );
		 // close file
		inFile.close();
	}
}


//////////////////////////////////////////////////////////////////////////

void
writeMoleculesToStream(	const sgm::Graph_Interface & multiMol
						, std::ostream & out
						, const bool ignoreProtons
						, const ggl::chem::GroupMap & groups )
{
	sgm::Graph_Interface::CompLabel compLabel;
	const size_t numOfEductMols = multiMol.connectedComponents(multiMol,compLabel);
	for (size_t i=0; i<numOfEductMols; ++i) {
		  // collect all nodes part of this molecule
		sgm::SubGraph::NodeList curNodes(compLabel.size());
		size_t numOfNodes = 0;
		for (size_t n=0; n<compLabel.size(); ++n) {
			  // check if this node is part of the molecule
			if (compLabel.at(n)==i) {
				curNodes[numOfNodes] = n;
				++numOfNodes;
			}
		}
		curNodes.resize(numOfNodes);
		  // get graph representation of this molecule
		sgm::SubGraph curSubGraph( multiMol, curNodes );
		  // copy to molecule instance
		Molecule curMol;
		MoleculeUtil::copy( curSubGraph, curMol );

		  // update SMIRKS
		if (i>0) { // add molecule separator if needed
			out <<".";
		}

		  // add SMILES of current molecule
		out << SMILESwriter::getSMILES(curMol, groups, ignoreProtons);
	}
}


//////////////////////////////////////////////////////////////////////////

