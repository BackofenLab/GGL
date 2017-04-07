
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

#include "ggl/Rule.hh"
#include "ggl/RuleGraph.hh"
#include "ggl/MR_ApplyRule.hh"

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
  /*! The used container to store the LeftSidePattern of Rules grouped by the
   * number of components of the LeftSidePattern.
   */
typedef std::vector< const sgm::Pattern_Interface* > RulePatternMap;

  //////////////////////////////////////////////////////////////////////////

typedef ggl::Graph GRAPHTYPE;

//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );  

//////////////////////////////////////////////////////////////////////////

void
parseRules(	std::istream & in
			, std::vector<ggl::Rule> & toFill
			, const size_t linesRead ) throw(std::exception);  

//////////////////////////////////////////////////////////////////////////

void
checkRules(	std::vector<ggl::Rule> & rules ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

void
parseGrid(	std::istream & in
			, GRAPHTYPE & toFill
			, size_t & numRow
			, size_t & numCol ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

void
printRule( std::ostream& out, ggl::Rule & rule );

//////////////////////////////////////////////////////////////////////////

void
printRules( std::ostream& out, std::vector<ggl::Rule> & rules );

//////////////////////////////////////////////////////////////////////////

void
printGrid(	std::ostream & out
			, const GRAPHTYPE & g
			, const size_t numRow
			, const size_t numCol );

//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample();

//////////////////////////////////////////////////////////////////////////

void
giveGridExample();

//////////////////////////////////////////////////////////////////////////

class MR_CopyMatch : public sgm::Match_Reporter
{

protected:
	  //! the graph to copy all matches to
	GRAPHTYPE & toFill;

	  //! access to the node labels of toFill
	boost::property_map< GRAPHTYPE, ggl::PropNodeLabel >::type nodeLabel;

	  //! number of rule applications done
	size_t changeCount;

public:

	MR_CopyMatch( GRAPHTYPE & toFill )
	 :	toFill(toFill),
	  	nodeLabel(boost::get( ggl::PropNodeLabel(), toFill )),
	  	changeCount(0)
	{}

	  //! forwards the match to the rule application handler and copies all
	  //! resulting changes into the toFill graph
	void
	reportHit (	const sgm::Pattern_Interface & pattern,
				const sgm::Graph_Interface & target,
				const sgm::Match & match )
	{

		// check for cast success
		assert( dynamic_cast<const ggl::LeftSidePattern*>(&pattern) != NULL /* pattern is no LeftSidePattern object ! */);
		// access to the rule that was matched
		const ggl::Rule & rule = static_cast<const ggl::LeftSidePattern*>(&pattern)->getRule();

		boost::property_map<	ggl::Rule::CoreGraph
								, ggl::Rule::NodeRightLabelProperty
					>::const_type ruleNodeRightLabel
						= boost::get(	ggl::Rule::NodeRightLabelProperty()
										, rule.getCore() );

		for (size_t n=0; n<match.size(); ++n) {
			nodeLabel[ boost::vertex( match.at(n), toFill ) ] = ruleNodeRightLabel[ boost::vertex( n, rule.getCore() ) ];
		}

		  // memorize rule application
		changeCount++;
	}

	  //! Access to the number of rule applications performed
	  //! @return the number of rule applications done
	size_t
	getChangeCount() const {
		return changeCount;
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
	InfoMode infoMode = OUT_NORMAL;
	
	std::istream* inRules = NULL;
	std::ifstream* inRulesFile = NULL;
	size_t inRulesLine = 0;
	std::istream* inGrid = &std::cin;
	std::ifstream* inGridFile = NULL;
	
	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;
	
	size_t iterations = 0;
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
	if (opts.argExist("ruleExample")) {
		giveRuleGMLExample();
		return 0;
	}
	if (opts.argExist("gridExample")) {
		giveGridExample();
		return 0;
	}
	if (opts.argExist("version")) {
		giveVersion();
		return 0;
	}
	if (opts.noErrors()) {
		  // rules are officially not mandatory.. we have to request manually
		if (!opts.argExist("rules")) {
			std::cerr <<"\n\n\tERROR : needed argument not given : 'rules' check usage\n\n";
			return -1;
		}
		if (opts.getBoolVal("v")) {
			infoMode = OUT_VERBOSE;
		}
	} else {
		return -1;
	}
	
	try {

		//////////////////////////////////////////////////////////////
		// set up iteration parameters
		//////////////////////////////////////////////////////////////

		if (opts.getIntVal("iter") <= 0) {
			std::ostringstream oss;
			oss	<<"number of rule application iterations ("
			<<opts.getIntVal("iter") <<") has to be at least 1";
			throw ArgException(oss.str());
		}
		iterations = (size_t)opts.getIntVal("iter");


		//////////////////////////////////////////////////////////////
		// set up input and output streams
		//////////////////////////////////////////////////////////////

		  // set Rule input stream
		if (opts.getStrVal("rules").size() == 0) {
			throw ArgException("no input file given");
		} else if (opts.getStrVal("rules").compare("STDIN") != 0) {
			inRulesFile = new std::ifstream( opts.getStrVal("rules").c_str()
										, std::ifstream::in );
			if (!inRulesFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open graph grammar rule input file '" <<opts.getStrVal("rules") <<"'";
				throw ArgException(oss.str());
			}
			inRules = inRulesFile;
		} else {
			  // read rules from STDIN
			inRules = &std::cin;
		}
		  // set GRID input stream
		if (opts.getStrVal("grid").size() == 0) {
			throw ArgException("no input grid file given");
		} else if (opts.getStrVal("grid").compare("STDIN") != 0) {
			inGridFile = new std::ifstream(	opts.getStrVal("grid").c_str()
										, std::ifstream::in );
			if (!inGridFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open input grid file '" <<opts.getStrVal("grid") <<"'";
				throw ArgException(oss.str());
			}
			inGrid = inGridFile;
		} else {
			  // read graphs from STDIN
			inGrid = &std::cin;
			  // check if rules already set to STDIN
			if (inGrid == inRules) {
				std::ostringstream oss;
				oss	<<"cannot read both, RULES and GRID, from STDIN";
				throw ArgException(oss.str());
			}
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
		parseRules( *inRules, rules, inRulesLine );

		checkRules( rules );

		  // generate left side pattern of each rule needed for its application
		  // and store separate rule lists based on the component number of
		  // their left side pattern
		for (size_t r=0; r<rules.size(); ++r) {
			ggl::LeftSidePattern* pattern = new ggl::LeftSidePattern(rules[r]);
			rulePattern.push_back( pattern );
		}
//		sgm::Graph_Interface::CompLabel compLabel;
//		for (size_t r=0; r<rules.size(); ++r) {
//			ggl::LeftSidePattern* pattern = new ggl::LeftSidePattern(rules[r]);
//			size_t compNumber = sgm::Graph_Interface::connectedComponents( *pattern, compLabel );
//			if (rulePattern.size() <= compNumber) {
//				rulePattern.resize(compNumber+1);
//			}
//			rulePattern[compNumber].push_back(pattern);
//		}

		GRAPHTYPE curGrid;
		GRAPHTYPE nextGrid;

		size_t numRow = 0, numCol = 0;

		  // parse the input grid
		parseGrid( *inGrid, curGrid, numRow, numCol );
		  // copy the input grid
		nextGrid = curGrid;


		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////

		  // set up graph matcher
		sgm::SGM_vf2 sgm;

		  // print initial grid
		(*out) <<"\n";
		printGrid( *out, curGrid, numRow, numCol );
		(*out) <<std::endl;

		  // perform iterations
		for (size_t it=0; it<iterations; ++it ) {

			  // set up Rule applyer
			MR_CopyMatch mr_applyRule( nextGrid );
			  // create target graph
			const sgm::Graph_boost< GRAPHTYPE > target( curGrid );

			  // FOR EACH set of rules with equal number of connected components
			  //          in leftsidepattern :
			  // find matches and apply rules
			  // apply rules
			sgm.findMatches( rulePattern
							, target
							, mr_applyRule
							, UINT_MAX );

			if (infoMode >= OUT_NORMAL) {
				  // output iteration count
				(*out) <<" " <<(it+1)<<". iteration\n" <<std::endl;
				  // print output grid
				printGrid( *out, nextGrid, numRow, numCol );
				(*out)	<<std::endl;
			}

			if (mr_applyRule.getChangeCount() == 0) {
				if (infoMode >= OUT_NORMAL) {
					(*out) <<"\n stopping iteration since no more changes possible ...\n" <<std::endl;
				}
				break;
			}

			  // set next grid to use
			curGrid = nextGrid;

		}

	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	inRules = &std::cin;
	inGrid = &std::cin;
	out = &std::cout;
	if (inRulesFile != NULL)		delete inRulesFile;
	if (inGridFile != NULL)		delete inGridFile;
	if (outFile != NULL)			delete outFile;
	for (size_t p=0; p<rulePattern.size(); ++p)
		delete rulePattern[p];

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
		"Explores the Game of Life for a given set of rules and a given"
		" start grid.\n"
		"\n"
		"The rule input has to be in GML format. The start grid has to be"
		" defined by a tabular format with label encoding."
		" Use '-ruleExample' or '-gridExample' for an example.\n"
		;
	
	allowedArgs.push_back(biu::COption(	
							"rules", false, biu::COption::STRING,
							"Graph grammar rules input file name or 'STDIN' when to read from standard input (use -ruleExample for a sample)"));
	allowedArgs.push_back(biu::COption(	
							"grid", false, biu::COption::STRING,
							"The initial grid input file name or 'STDIN' when to read from standard input (see -gridExample for a sample)",
							"STDIN"));
	allowedArgs.push_back(biu::COption(	
							"iter", true, biu::COption::INT, 
							"Number of rule application iterations",
							"1"));
	allowedArgs.push_back(biu::COption(	
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(	
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"ruleExample", true, biu::COption::BOOL, 
							"Displays an example for the graph grammar rule GML encoding"));
	allowedArgs.push_back(biu::COption(	
							"gridExample", true, biu::COption::BOOL,
							"Displays an example for the input grid encoding"));
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
 * Parses the given input stream for rules and adds each rule to the given
 * container.
 * 
 * Each rule has to be headed by a comment line starting with a '#' character 
 * that specifies the name of the rule.
 * 
 * @param in the input stream to read from
 * @param toFill the container to add the rules to
 * @param linesRead the number of lines already read from input (needed for 
 *                  error reporting)
 * 
 */
void
parseRules(	std::istream & in
			, std::vector<ggl::Rule> & toFill
			, const size_t linesRead = 0 )  throw(std::exception)
{
	std::string line;
	size_t lineNumber = linesRead;
	std::string ruleString = std::string("");
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;
		
		  // check if next rule header or end of file was reached
		if (line[0] == '#' || in.eof()) {
			  // convert rule string to rule object
			if (ruleString.size() > 0) {
				std::pair<ggl::Rule, int> 
						ret = ggl::Rule_GMLparser::parseRule( ruleString );
				  // check for parsing error
				if (ret.second != -1) {
					std::ostringstream oss;
					oss	<<" Rule parsing error in rule specification BEFORE line "
						<<line
						<<" at rule string position " <<ret.second;
					throw ArgException(oss.str());
				}
				toFill.push_back(ret.first);
			}
			  // clear processed rule
			ruleString.clear();
		} else {
			  // append line content to rule string
			ruleString.append(line);
		}
	}
}


//////////////////////////////////////////////////////////////////////////


void
checkRules(	std::vector<ggl::Rule> & rules ) throw(std::exception)
{
	for (size_t r=0; r<rules.size(); ++r) {
		if ( boost::num_vertices( rules.at(r).getCore() ) != 1 ) {
			std::stringstream err;
			err <<"ERROR : rule " <<rules.at(r).getID() <<" : "
					<<"rule size has to be one!\n";
			throw ArgException(err.str());
		}
		boost::property_map< ggl::Rule::CoreGraph
								, ggl::Rule::NodeContextProperty
							>
			::const_type ruleNodeContext = boost::get(
					ggl::Rule::NodeContextProperty()
								, rules.at(r).getCore() );
		if ( ruleNodeContext[ boost::vertex( 0, rules.at(r).getCore() ) ] != ggl::Rule::RULE_LABEL_CHANGE ) {
			std::stringstream err;
			err <<"ERROR : rule " <<rules.at(r).getID() <<" : "
					<<"only label change allowed!\n";
			throw ArgException(err.str());
		}
	}
}


//////////////////////////////////////////////////////////////////////////

void
parseGrid(	std::istream & in
			, GRAPHTYPE & toFill
			, size_t & numRow
			, size_t & numCol ) throw(std::exception)
{
	size_t lineCount = 0;

	if (boost::num_vertices( toFill ) > 0) {
		throw ArgException("input grid graph is not empty");
	}

	numRow = 0;
	numCol = 0;

	std::vector< std::vector< int > > labels;
	std::map< int, std::string > int2string;

	while (in.good()) {
		  // read next line
		std::string line = "";
		std::getline( in, line );
		lineCount++;

		  // stop if empty line
		if (line.find_first_not_of(" \t\n") == std::string::npos) {
			break;
		}
		  // parse line
		std::stringstream ss(line);
		std::vector<int> row;
		int nextLabel= -1;
		  // read until no 0 or 1 can be parsed
		while( ss.good() ) {
			nextLabel = -1;
			ss >> nextLabel;
			if (nextLabel == 0 || nextLabel == 1) {
				row.push_back(nextLabel);
				  // store string label if unknown
				if (int2string.find(nextLabel)==int2string.end()) {
					std::stringstream conv;
					conv <<nextLabel;
					int2string[nextLabel] = conv.str();
				}
			} else {
				break;
			}
		}
		  // check if parse was ok
		if (labels.size() == 0) {
			if (row.size() == 0) {
				throw ArgException(" Grid parsing error in line 1 : no entries in first row");
			}
		} else
		if (row.size() != labels.rbegin()->size()) {
			std::ostringstream oss;
			oss	<<" Grid parsing error in line "
				<<lineCount
				<<" : entry number differs to preceding rows";
			throw ArgException(oss.str());
		}
		  // store parsed row labels
		labels.push_back(row);
		numRow++;

	}

	numCol = labels.begin()->size();

	if (numCol == 0 || numRow == 0) {
		throw ArgException("ERROR : no input grid could be parsed");
	}

	if (numCol < 3 || numRow < 3 ) {
		throw ArgException("ERROR : minimal size of input grid is 3x3");
	}


	  // create graph instance from parsed label matrix
	boost::property_map< GRAPHTYPE, ggl::PropNodeLabel >::type nodeLabel = boost::get( ggl::PropNodeLabel(), toFill );
	boost::property_map< GRAPHTYPE, ggl::PropEdgeLabel >::type edgeLabel = boost::get( ggl::PropEdgeLabel(), toFill );

	  // insert nodes
	GRAPHTYPE::vertex_descriptor newVertex;
	for (size_t r=0; r<numRow; ++r) {
		for (size_t c=0; c<numCol; ++c) {
			  // add new vertex
			newVertex = boost::add_vertex( toFill );
			nodeLabel[newVertex] = int2string[ labels[r][c] ];
		}
	}
	  // add row edges
	std::pair< GRAPHTYPE::edge_descriptor, bool > newEdge;
	for (size_t r=0; r<numRow; ++r) {
		size_t nodeShift = r*numCol;
		for (size_t c=0; c<numCol; ++c) {
			size_t neigh = (c==(numCol-1) ? 0 : c+1 );
			newEdge = boost::add_edge( boost::vertex( nodeShift+c, toFill )
									,  boost::vertex( nodeShift+neigh, toFill )
									, toFill );
			edgeLabel[ newEdge.first ] = "-";
		}
	}
	  // add col edges
	for (size_t c=0; c<numCol; ++c) {
		for (size_t r=0; r<numRow; ++r) {
			size_t neigh = (r==(numRow-1) ? 0 : r+1 );
			newEdge = boost::add_edge( boost::vertex( c+(r*numCol), toFill )
									,  boost::vertex( c+(neigh*numCol), toFill )
									, toFill );
			edgeLabel[ newEdge.first ] = "|";
		}
	}
	  // add diagonal edges
	for (size_t c=0; c<numCol; ++c) {
		for (size_t r=0; r<numRow; ++r) {

			size_t idx_cr = (r*numCol) + c;

			int rowNeigh = (r==(numRow-1) ? (-r*numCol) : +numCol );
			int colNeigh = (c==(numCol-1) ? (-c) : +1 );

			newEdge = boost::add_edge( boost::vertex( idx_cr, toFill )
									,  boost::vertex( idx_cr + rowNeigh + colNeigh, toFill )
									, toFill );
			edgeLabel[ newEdge.first ] = "\\";
			newEdge = boost::add_edge( boost::vertex( idx_cr + colNeigh, toFill )
									,  boost::vertex( idx_cr + rowNeigh, toFill )
									, toFill );
			edgeLabel[ newEdge.first ] = "/";
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
printGrid(	std::ostream & out
			, const GRAPHTYPE & g
			, const size_t numRow
			, const size_t numCol )
{
	if (boost::num_vertices( g ) == 0) {
		return;
	}

	boost::property_map< GRAPHTYPE, ggl::PropNodeLabel >::const_type nodeLabel = boost::get( ggl::PropNodeLabel(), g );


	for (size_t c=0; c<=(numCol+1); ++c)
		out <<"#";
	out <<"\n";
	for (size_t r=0; r<numRow; ++r) {
		size_t nodeShift = r*numCol;
		out <<"#";
		for( size_t c=0; c<numCol; ++c ) {
//			out << (c>0?" ":"") <<nodeLabel[ boost::vertex( nodeShift+c, g ) ];
			out << (nodeLabel[ boost::vertex( nodeShift+c, g ) ] == "0" ? " " : "O");
		}
		out <<"#\n";
	}
	for (size_t c=0; c<=(numCol+1); ++c)
		out <<"#";
	out <<"\n";
}


//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample()
{
	std::cout <<"\n"
"\n"
"rule [\n"
" ruleID \"3 neighbors = 1\"\n"
" left [\n"
"   node [ id 0 label \"0\" ]\n"
" ]\n"
" right [\n"
"   node [ id 0 label \"1\" ]\n"
" ]\n"
" constrainAdj [\n"
"   id 0\n"
"   op =\n"
"   count 3\n"
"   nodeLabels [\n"
"     label \"1\"\n"
"   ]\n"
" ]\n"
"]\n"
"\n"
	<<std::endl;
}

//////////////////////////////////////////////////////////////////////////

void
giveGridExample()
{
	std::cout <<""
"1 0 0 0 0 0 0 0\n"
"0 1 0 0 0 0 0 0\n"
"1 0 0 0 0 0 0 0\n"
"0 0 0 1 0 0 0 1\n"
"1 0 0 0 0 0 0 0\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////


