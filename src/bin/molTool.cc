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

#include <sgm/GraphScaffold.hh>

#include <ggl/Graph.hh>
#include <ggl/Graph_GML_writer.hh>
#include <ggl/Graph_gSpan_writer.hh>
#include <ggl/chem/Molecule.hh>
#include <ggl/chem/SMILESwriter.hh>

#include <ggl/chem/AP_disabled.hh>


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
  * pushes the scaffold annotation as class information into the atom labels
  */
void
addScaffoldAnnotation( ggl::chem::Molecule & mol );

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

	AromaticityPerception * aromaticityPrediction = NULL;

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
		if (boost::iequals(opts.getStrVal("in"),"STDIN")) {
			  // read from STDIN
			in = &std::cin;
		} else {
			  // read from file
			inFile = new std::ifstream( opts.getStrVal("in").c_str()
										, std::ifstream::in );
			if (!inFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open input file '" <<opts.getStrVal("in") <<"'";
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


		// check outMode parameter setup : non-empty and within allowed chars
		const std::string allowedModeChars = "sSgGlLpP";
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
		ggl::chem::GroupMap * outputGroups = &noGroups;
		if ( ! opts.argExist("noGroupInsertion")) {
			outputGroups = &moleculeGroups;
		}

		SMILES_container parsedMols;
		if (opts.getCharVal("inMode") == 'S' || opts.getCharVal("inMode") == 's') {
			parseSMILES( *in, parsedMols, inLine, moleculeGroups );
		} else if (opts.getCharVal("inMode") == 'G' || opts.getCharVal("inMode") == 'g') {
			parseMolGML( *in, parsedMols, inLine, moleculeGroups );
		} else {
			std::ostringstream oss;
			oss	<<"Input format '" <<opts.getStrVal("inMode") <<"' is unknown";
			throw ArgException(oss.str());
		}

		  // do molecule correction
		if (!opts.getBoolVal("noInputCorrection")) {

			  // set aromaticity model
			switch(opts.getCharVal("aromaticity")) {
			case 'o' :
			case 'O' :
				aromaticityPrediction = new AP_NSPDK("OpenBabel:2013");
				break;
			case 'm' :
			case 'M' :
				aromaticityPrediction = new AP_NSPDK("Marvin:general:2013");
				break;
			case 'n' :
			case 'N' :
				aromaticityPrediction = new AP_disabled();
				break;
			default:
					std::ostringstream oss;
					oss	<<"aromaticity perception type '" <<opts.getCharVal("aromaticity") <<"'"
						<<" is not supported";
					throw ArgException(oss.str());
			}

			 // ensure protons are represented if aromaticity is to be predicted
			if (opts.argExist("noProtonFilling") && aromaticityPrediction != NULL) {
				throw ArgException("Aromaticity perception requires proton filling!");
			}

			// correct aromaticity of input molecules
			SMILES_container producedSmiles;
			correctInputMolecules( parsedMols, producedSmiles, aromaticityPrediction,!opts.argExist("noProtonFilling") );
			  // clear input molecules
			for (SMILES_container::iterator it=parsedMols.begin(); it!=parsedMols.end() ; ++it ) {
				delete it->second;
			}
			parsedMols.clear();
			  // copy data back
			parsedMols = producedSmiles;
		}

		  // do molecule sanity checks
		if (!opts.getBoolVal("noInputCheck")) {
			for (SMILES_container::iterator it=parsedMols.begin(); it!=parsedMols.end() ; ++it ) {
				size_t conStatus = ggl::chem::MoleculeUtil::isConsistent( *(it->second) );
				if( conStatus != ggl::chem::MoleculeUtil::C_Consistent ) {
					(*out) <<"\n PROBLEM : molecule '"
							<<it->first <<"' is not chemically correct"
								" or contains unsupported properties:\n";
					ggl::chem::MoleculeUtil::decodeConsistencyStatus( conStatus, (*out) );
				}
			}
		}



		//////////////////////   OUTPUT  /////////////////////////////////////

		  // write results to stream
		if (opts.getCharVal("outMode") == 'S' || opts.getCharVal("outMode") == 's') {
			  // write in SMILES format
			for( SMILES_container::const_iterator m=parsedMols.begin();
					m != parsedMols.end(); ++m )
			{
				ggl::chem::Molecule tmp;
				ggl::chem::MoleculeUtil::copy(*(m->second), tmp);
				  // add annotation if necessary
				if (opts.getBoolVal("scaffoldAnnotation")) {
					addScaffoldAnnotation(tmp);
				}
				try {
					if (opts.getBoolVal("noProtonRemoval")) {
						  // get proton pruned molecule
						*out <<ggl::chem::SMILESwriter::getSMILES(tmp, *outputGroups,false) <<std::endl;
					} else {
						  // write molecule SMILES to stream, protons are put into complex atom labels
						*out <<ggl::chem::SMILESwriter::getSMILES(tmp, *outputGroups, true) <<std::endl;
					}
				} catch (std::exception & ex) {
					std::cerr <<"ERROR in SMILES generation for molecule '" <<m->first <<"' : " <<ex.what() <<std::endl;
				}

			}
		} else if (opts.getCharVal("outMode") == 'G' || opts.getCharVal("outMode") == 'g'
				|| opts.getCharVal("outMode") == 'L' || opts.getCharVal("outMode") == 'l'
				|| opts.getCharVal("outMode") == 'P' || opts.getCharVal("outMode") == 'p')
		{
			  // write in GML format
			for( SMILES_container::const_iterator m=parsedMols.begin();
					m != parsedMols.end(); ++m )
			{
				ggl::chem::Molecule tmp;
				ggl::chem::MoleculeUtil::copy(*(m->second), tmp);

				  // compress groups if available
				if (!outputGroups->empty()) {
					ggl::chem::MoleculeUtil::compressGroups(tmp, *outputGroups);
				}
				  // get proton pruned molecule if needed
				if (!opts.getBoolVal("noProtonRemoval")) {
					ggl::chem::MoleculeUtil::removeProtons(tmp);
				}
				  // add annotation if necessary
				if (opts.getBoolVal("scaffoldAnnotation")) {
					addScaffoldAnnotation(tmp);
				}

				  // create OUTPUT
				if (opts.getCharVal("outMode") == 'G' || opts.getCharVal("outMode") == 'g') {
					  // write molecule in GML to stream
					ggl::Graph_GML_writer::write( *out, tmp, true );
				} else if (opts.getCharVal("outMode") == 'L' || opts.getCharVal("outMode") == 'l') {
					  // write molecule in one line GML to stream
					ggl::Graph_GML_writer::write( *out, tmp, false );
				} else {
					  // write molecule in gSpan to stream
					ggl::Graph_gSpan_writer::write( *out, tmp );
				}
				*out <<std::endl;
			}
		} else {
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

	in = &std::cin;
	out = &std::cout;
	if (inFile != NULL)	{ inFile->close(); delete inFile; }
	if (outFile != NULL)		{ outFile->close(); delete outFile; }
	if (aromaticityPrediction != NULL)	delete aromaticityPrediction;

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
		"Enables the conversion of molecule graph into different formats link SMILES, GML, or gSpan.\n"
		" Note, a UNIQUE list of molecules is returned in case several input molecules encode the same compound.\n"
		" Furthermore, it does sanity checks and molecule corrections (as aromaticity perception, proton filling, etc.) for each input molecule.\n"
		" On top, a scaffold annotation can be added to the molecules.\n"
		;

	allowedArgs.push_back(biu::COption(
							"in", true, biu::COption::STRING,
							"Input file name or 'STDIN' when to read from standard input",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"inMode", true, biu::COption::CHAR,
							"Input format : (S)miles or (G)ML",
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
							"Output format : (S)miles, (G)ML, GML one (L)ine, or gS(p)an",
							"S"));
	allowedArgs.push_back(biu::COption(
							"scaffoldAnnotation", true, biu::COption::BOOL,
							"Add scaffold annotation (ring = 1, linker = 2, side chain = 3) to atom labels as class information."));
	allowedArgs.push_back(biu::COption(
							"noInputCorrection", true, biu::COption::BOOL,
							"Dont correct the input molecules (aromaticity perception, proton filling, ...)"));
	allowedArgs.push_back(biu::COption(
							"noProtonFilling", true, biu::COption::BOOL,
							"Dont add protons to the molecule before correcting it"));
	allowedArgs.push_back(biu::COption(
							"noProtonRemoval", true, biu::COption::BOOL,
							"Dont remove protons before producing the output"));
	allowedArgs.push_back(biu::COption(
							"noGroupInsertion", true, biu::COption::BOOL,
							"Dont replace molecule groups with their identifiers (from 'groups' parameter) before producing the output"));
	allowedArgs.push_back(biu::COption(
							"noInputCheck", true, biu::COption::BOOL,
							"Dont check the input molecules for consistency (atom/bond label, ...)"));
	allowedArgs.push_back(biu::COption(
							"aromaticity", true, biu::COption::CHAR,
							"The aromaticity perception model to be used : (M)arvin general model, (O)penBabel model, or (N)o aromaticity perception.",
							"M"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////

void
addScaffoldAnnotation( ggl::chem::Molecule & mol )
{
	  // create scaffold annotation object
	sgm::GraphScaffold scaffold;

	const sgm::GraphScaffold::ScaffoldAnnotation annotation
		= scaffold.getScaffoldAnnotation( ggl::chem::Molecule_Graph(mol) );

	ggl::chem::Molecule::vertex_iterator vIt, vItEnd;
	  // label access
	boost::property_map< ggl::chem::Molecule, ggl::chem::PropNodeLabel >::type
		nodeLabel = boost::get( ggl::chem::PropNodeLabel(), mol );
	  // index access
	boost::property_map< ggl::chem::Molecule, ggl::chem::PropNodeIndex >::const_type
		nodeIndex = boost::get( ggl::chem::PropNodeIndex(), mol );

	  // update all vertices
	for (boost::tie(vIt,vItEnd) = boost::vertices(mol); vIt != vItEnd; ++vIt) {
		  // get current atom label
		const std::string atomLabel = nodeLabel[*vIt];

		const std::string atom = ggl::chem::MoleculeUtil::getAtom( atomLabel );
		const size_t protonNumber = ggl::chem::MoleculeUtil::getProtons( atomLabel );
		const int charge = ggl::chem::MoleculeUtil::getCharge( atomLabel );

		size_t classID = 0;
		switch ( annotation.at(nodeIndex[*vIt])) {
			case sgm::GraphScaffold::GST_RING:
				classID = 1;
				break;
			case sgm::GraphScaffold::GST_LINKER:
				classID = 2;
				break;
			case sgm::GraphScaffold::GST_DANGLING:
				classID = 3;
				break;
			default:
				throw ArgException("found node with unknown scaffold annotation");
				break;
		}

		nodeLabel[ *vIt ] = ggl::chem::MoleculeUtil::getComplexAtomLabel( atom, protonNumber, charge, classID );
	}


}


//////////////////////////////////////////////////////////////////////////
