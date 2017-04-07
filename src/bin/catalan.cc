
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "biu/OptionParser.hh"

#include "sgm/GM_vf2.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"
#include "sgm/MR_stream.hh"
#include "sgm/MR_SymmBreak.hh"
#include "sgm/PA_OrderCheck.hh"
#include "sgm/Graph_boost.hh"


#include "ggl/Rule.hh"
#include "ggl/RuleGraph.hh"
#include "ggl/MR_ApplyRule.hh"
#include "ggl/DFS_ApplyRule.hh"
#include "ggl/GS_STL.hh"
#include "ggl/Graph_GMLparser.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/GS_stream.hh"

#include "version.hh"


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


  //////////////////////////////////////////////////////////////////////////
  /*! The used container to store the LeftSidePattern of Rules
   */
typedef std::vector< const sgm::Pattern_Interface* > RulePatternMap;


  //////////////////////////////////////////////////////////////////////////

typedef ggl::Graph GRAPHTYPE;

typedef std::map< std::string, GRAPHTYPE > Graph_container;
typedef std::vector< const GRAPHTYPE * > GraphP_container;

//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );  

//////////////////////////////////////////////////////////////////////////

void
parseRules(	std::vector<ggl::Rule> & toFill );

//////////////////////////////////////////////////////////////////////////

void
parseGraph(	std::istream & in
			, GRAPHTYPE & toFill
			, const size_t linesRead ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

void
printRule( std::ostream& out, ggl::Rule & rule );

//////////////////////////////////////////////////////////////////////////

void
printRules( std::ostream& out, std::vector<ggl::Rule> & rules );


//////////////////////////////////////////////////////////////////////////

void
giveInputExample();

//////////////////////////////////////////////////////////////////////////

typedef ggl::DFS_ApplyRule Catalan_DFS;

class CatalanVisitor : public Catalan_DFS::DFS_Visitor {

public:

	virtual
	Decision
	status( const GRAPHTYPE & graph ) {
		  // check if a solution was found
		  // i.e. only a single node is left
		if (boost::num_vertices(graph) == 1 && boost::num_edges(graph) == 0) {
			std::cout <<"\n Found a solution\n\n" <<std::endl;
			  // only one solution wanted
			return SOLUTION_STOP;
		}
//// DEBUG OUTPUT
//ggl::Graph_GML_writer gml_writer;
//gml_writer.write( std::cout, graph, true );

		  // standard case : continue search
		return CONTINUE;
	}

};



//////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv ) {
	using namespace std;
	
	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////
	int exitValue = 0;
	
	enum InfoMode {OUT_SILENT, OUT_NORMAL, OUT_VERBOSE};
//	InfoMode infoMode = OUT_NORMAL;
	
	std::istream* inGraph = &std::cin;
	std::ifstream* inGraphFile = NULL;
	
	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;
	
	 // rule pattern for each number of connected components
	RulePatternMap rulePattern;
	
	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////
	
	biu::OptionMap allowedArgs;	//< map of allowed arguments
	std::string infoText;		//< info string of the program
	initAllowedArguments(allowedArgs,infoText);	// init
	
		// parse programm arguments	
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
													argv, infoText);
		// arguments parseable and all mandatory arguments given
		// help output
	if (opts.argExist("help")) {
		opts.coutUsage();
		return 0;
	}
	if (opts.argExist("inputExample")) {
		giveInputExample();
		return 0;
	}
	if (opts.argExist("version")) {
		giveVersion();
		return 0;
	}
	if (opts.noErrors()) {
//		if (opts.getBoolVal("v")) {
//			infoMode = OUT_VERBOSE;
//		}
	} else {
		return -1;
	}
	
	try {


		//////////////////////////////////////////////////////////////
		// set up input and output streams
		//////////////////////////////////////////////////////////////

		  // set GRID input stream
		if (opts.getStrVal("catalan").size() == 0) {
			throw ArgException("no input grid file given");
		} else if (opts.getStrVal("catalan").compare("STDIN") != 0) {
			inGraphFile = new std::ifstream(	opts.getStrVal("catalan").c_str()
										, std::ifstream::in );
			if (!inGraphFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open input grid file '" <<opts.getStrVal("catalan") <<"'";
				throw ArgException(oss.str());
			}
			inGraph = inGraphFile;
		} else {
			  // read graphs from STDIN
			inGraph = &std::cin;
		}
		  // set output stream
		if (opts.getStrVal("out").size() == 0) {
			throw ArgException("no output file given");
		} else if (opts.getStrVal("out").compare("STDOUT") != 0) {
			outFile = new std::ofstream(	opts.getStrVal("out").c_str()
											, std::ofstream::out );
			if (!outFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open output file '" <<opts.getStrVal("out") <<"'";
				throw ArgException(oss.str());
			}
			out = outFile;
		}


		//////////////////////////////////////////////////////////////
		// parse Rules and Graphs from input
		//////////////////////////////////////////////////////////////

		std::vector<ggl::Rule> rules;
		parseRules( rules );

//		printRules(*out, rules);

		  // generate left side pattern of each rule needed for its application
		  // and store separate rule lists based on the component number of
		  // their left side pattern
		for (size_t r=0; r<rules.size(); ++r) {
			ggl::LeftSidePattern* pattern = new ggl::LeftSidePattern(rules[r]);
			rulePattern.push_back( pattern );
		}


		  // parse the input grid
		GRAPHTYPE inputCatalan;
		size_t linesRead = 0;
		parseGraph( *inGraph, inputCatalan, linesRead );

		*out <<"\n Input Catalan \n\n"
			<<sgm::Graph_boost<GRAPHTYPE>(inputCatalan)
			<<std::endl;

		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////


		  // set up graph matcher
		sgm::SGM_vf2 sgm;

		Catalan_DFS::SearchTrace dfs_trace;

		Catalan_DFS catalan_dfs;
		CatalanVisitor visitor;

		 // setup dummy graph storage for the solution
		class DoNothing : public ggl::Graph_Storage {
		public: virtual void add( const GRAPHTYPE & graph ) {};
		} solutionStorage;

		  // perform DFS
		bool solutionFound = catalan_dfs.findSolution( rules, inputCatalan, solutionStorage, visitor, &dfs_trace, true, &sgm );


		if ( !solutionFound ) {
			*out <<" no solution found ... \n"<<std::endl;
		} else {
			ggl::Graph_GML_writer gml_writer;
			  // print solution trace
			*out <<" --> the DFS trace :\n" <<std::endl;;
			for (size_t i=0; i<dfs_trace.size(); ++i) {
				*out <<"\n";
				gml_writer.write( *out, *(dfs_trace.at(i)), true );
				*out <<std::endl;
			}
		}

		for (size_t i=0; i<dfs_trace.size(); ++i) {
			delete dfs_trace[i];
		}
		dfs_trace.clear();


	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	inGraph = &std::cin;
	out = &std::cout;
	if (inGraphFile != NULL)		delete inGraphFile;
	if (outFile != NULL)			delete outFile;
	for (size_t i=0; i<rulePattern.size(); ++i ) {
		delete rulePattern[i];
	}

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
		"Solves a Catalan via depths-first-search given a start setup.\n"
		"\n"
		"The start setup has to be defined by a tabular format."
		" Use '-inputExample' for an example.\n"
		;
	
	allowedArgs.push_back(biu::COption(	
							"catalan", false, biu::COption::STRING,
							"The initial catalan input file name or 'STDIN' when to read from standard input (see -inputExample for a sample)",
							"STDIN"));
	allowedArgs.push_back(biu::COption(	
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"inputExample", true, biu::COption::BOOL,
							"Displays an example for the Catalan input encoding"));
	allowedArgs.push_back(biu::COption(	
							"v", true, biu::COption::BOOL, 
							"Verbose output"));
	allowedArgs.push_back(biu::COption(	
							"version", true, biu::COption::BOOL, 
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////

#include "ggl/Rule_GMLparser.hh"

/*!
 * Create the catalan rules
 * 
 * @param toFill the container to add the rules to
 * 
 */
void
parseRules(	std::vector<ggl::Rule> & toFill )
{
	{
		std::string pruneRuleString =
				"rule [ \n"
				" ruleID \"loop cleanup\" \n"
				" wildcard \"*\" \n"
				" context [ \n"
				"   node [ id 0 label \"*\" ] \n"
				" ] \n"
				" left [ \n"
				"   edge [ source 0 target 0 label \"*\" ] \n"
				" ] \n"
				"]";
		  // parse rule
		std::pair<ggl::Rule, int>
			ret = ggl::Rule_GMLparser::parseRule( pruneRuleString );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in prune rule specification "
				<<" at rule string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store rule
		toFill.push_back(ret.first);
	}
	{
		std::string pruneRuleString =
				"rule [ \n"
				" ruleID \"parallel edge cleanup\" \n"
				" wildcard \"*\" \n"
				" context [ \n"
				"   node [ id 0 label \"*\" ] \n"
				"   node [ id 1 label \"*\" ] \n"
				"   edge [ source 0 target 1 label \"*\" ] \n"
				" ] \n"
				" left [ \n"
				"   edge [ source 0 target 1 label \"*\" ] \n"
				" ] \n"
				" constrainNoEdge [ source 0 target 0 ]\n"
				" constrainNoEdge [ source 1 target 1 ]\n"
				"]";
		  // parse rule
		std::pair<ggl::Rule, int>
			ret = ggl::Rule_GMLparser::parseRule( pruneRuleString );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in prune rule specification "
				<<" at rule string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store rule
		toFill.push_back(ret.first);
	}
	{
		std::string shrinkRuleString =
					"rule [ \n"
					" ruleID \"catalan shrink\" \n"
					" wildcard \"*\" \n"
					" left [ \n"
					"   node [ id 0 label \"*\" ] \n"
					"   node [ id 1 label \"*\" ] \n"
					"   node [ id 2 label \"*\" ] \n"
					"   node [ id 3 label \"*\" ] \n"
					"   edge [ source 0 target 1 label \"*\" ] \n"
					"   edge [ source 0 target 2 label \"*\" ] \n"
					"   edge [ source 0 target 3 label \"*\" ] \n"
					" ] \n"
					" right [ \n"
					"   node [ id 4 label \"N\" ] \n"
					" ] \n"
					" constrainNoEdge [ source 0 target 0 ]\n"
					" constrainAdj [ \n"
					"   id 0 \n"
					"   op = \n"
					"   count 3 \n"
					"   edgeLabels [ label \"*\" ] \n"
					"   nodeLabels [ label \"*\" ] \n"
					" ] \n"
					" copyAndPaste [ source 1 id 4 edgeLabels [ label \"*\" ] ] \n"
					" copyAndPaste [ source 2 id 4 edgeLabels [ label \"*\" ] ] \n"
					" copyAndPaste [ source 3 id 4 edgeLabels [ label \"*\" ] ] \n"
					"]";

		  // parse rule
		std::pair<ggl::Rule, int>
				ret = ggl::Rule_GMLparser::parseRule( shrinkRuleString );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in shrink rule specification"
				<<" at rule string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store rule
		toFill.push_back(ret.first);
	}
}


//////////////////////////////////////////////////////////////////////////

void
parseGraph(	std::istream & in
			, GRAPHTYPE & toFill
			, const size_t linesRead ) throw(std::exception)
{
	std::string line;
	size_t lineNumber = linesRead;
	std::string graphString = std::string("");
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;

		  // check if next rule header or end of file was reached
		if (line[0] == '#' || in.eof()) {
			  // convert rule string to rule object
			if (graphString.size() > 0) {
				std::pair<ggl::Graph, int>
						ret = ggl::Graph_GMLparser::parseGraph( graphString );
				  // check for parsing error
				if (ret.second != -1) {
					std::ostringstream oss;
					oss	<<" Graph parsing error in graph specification BEFORE line "
						<<line
						<<" at graph string position " <<ret.second;
					throw ArgException(oss.str());
				}
				toFill = ret.first;
				return;
			}
			  // clear processed rule
			graphString.clear();
		} else {
			  // append line content to rule string
			graphString.append(line);
		}
	}
}


//////////////////////////////////////////////////////////////////////////


void
printRules( std::ostream& out, std::vector<ggl::Rule> & rules )
{
	for (size_t i=0; i<rules.size(); ++i ) {
		printRule(out, rules[i]);
		out <<"\n";
	}

}


//////////////////////////////////////////////////////////////////////////


void
printRule( std::ostream& out, ggl::Rule & rule )
{
	{
		out <<"\n == LEFT_SIDE_PATTERN ==\n" <<std::endl;
		ggl::LeftSidePattern g(rule);
		for (size_t i=0; i<g.getNodeNumber(); ++i ) {
			out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
					eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
			{
				out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			out <<" |\n";
		}
	}
	{
		out <<"\n == RIGHT_SIDE_PATTERN ==\n" <<std::endl;
		ggl::RightSidePattern g(rule);
		for (size_t i=0; i<g.getNodeNumber(); ++i ) {
			out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
					eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
			{
				out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			out <<" |\n";
		}
	}
	out <<"\n";
}



//////////////////////////////////////////////////////////////////////////

void
giveInputExample()
{
	std::cout <<""
"graph [ \n"
" node [ id 1 label \"1\" ]\n"
" node [ id 2 label \"2\" ]\n"
" node [ id 3 label \"3\" ]\n"
" node [ id 4 label \"4\" ]\n"
" node [ id 5 label \"5\" ]\n"
" node [ id 6 label \"6\" ]\n"
" node [ id 7 label \"7\" ]\n"
" edge [ source 1 target 2 label \"1-2\" ]\n"
" edge [ source 2 target 3 label \"2-3\" ]\n"
" edge [ source 2 target 4 label \"2-4\" ]\n"
" edge [ source 3 target 5 label \"1-2\" ]\n"
" edge [ source 4 target 5 label \"1-2\" ]\n"
" edge [ source 5 target 6 label \"1-2\" ]\n"
" edge [ source 5 target 7 label \"1-2\" ]\n"
"]\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////


