
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

#include "ggl/Rule.hh"
#include "ggl/RuleGraph.hh"
#include "ggl/MR_ApplyRule.hh"
#include "ggl/DFS_ApplyRule.hh"
#include "ggl/GS_STL.hh"
//#include "ggl/GS_stream.hh"

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
giveInputExample();

//////////////////////////////////////////////////////////////////////////

class GS_GridString : public ggl::Graph_Storage {
	Graph_container & storage;
public:
	GS_GridString( Graph_container& storage )
	 :	storage(storage)
	{}

	void
	add( const GRAPHTYPE & graph )
	{
		  // get grid string representation
		std::string str = getGridString( graph );
		  // store if not already present
		if (storage.find(str) == storage.end()) {
			storage[str] = graph;
		}
	}

	static
	std::string
	getGridString( const GRAPHTYPE& graph ) {
		  // get grid node labels in order
		std::vector< std::string > grid = getGridOrdered( graph );
		  // generate string representation of grid
		std::stringstream s;
		for (size_t i=0; i<grid.size(); ++i ) {
			s <<grid.at(i) <<":";
		}
		  // return string representation
		return s.str();
	}

	static
	std::vector< std::string >
	getGridOrdered( const GRAPHTYPE& g )
	{
		std::vector< std::string > labelInOrder(boost::num_vertices(g),"");

		boost::property_map< GRAPHTYPE, ggl::PropEdgeLabel >::const_type edgeLabel = boost::get( ggl::PropEdgeLabel(), g );
		boost::property_map< GRAPHTYPE, ggl::PropNodeLabel >::const_type nodeLabel = boost::get( ggl::PropNodeLabel(), g );

		GRAPHTYPE::vertex_descriptor v;
		GRAPHTYPE::out_edge_iterator a, a_end;

		for (size_t i=0; i<labelInOrder.size(); ++i ) {

			v = boost::vertex( i, g );
			size_t count = 0;
			int min = (int)labelInOrder.size();
			int max = -1;
			for (boost::tie(a,a_end) = boost::out_edges( v, g ); a != a_end; ++a) {
				if (edgeLabel[ *a ] != "-") {
					count++;
					int tmp = (int)(atoi( edgeLabel[ *a ].c_str() ));
					if (tmp > max) { max = tmp; }
					if (tmp < min) { min = tmp; }
				}
			}
			  // handle last cell
			if (max == (int)(labelInOrder.size()-2) && min < (int)(labelInOrder.size()-3) ) {
				max = (int)(labelInOrder.size()-1);
			}
			  // write at correct position
			labelInOrder[max] = nodeLabel[ v ];

		}
		return labelInOrder;
	}

};


//////////////////////////////////////////////////////////////////////////

typedef ggl::DFS_ApplyRule Sudoku_DFS;

class SudokuVisitor : public Sudoku_DFS::DFS_Visitor {

public:

	virtual
	Decision
	status( const GRAPHTYPE & graph ) {
		  // check if graph is completely filled and a solution was found
		bool allFilled = true;
		boost::property_map< GRAPHTYPE, ggl::PropNodeLabel >::const_type nodeLabel = boost::get( ggl::PropNodeLabel(), graph );
		for (int i=boost::num_vertices(graph)-1; allFilled && i>=0; i--) {
			allFilled = nodeLabel[ boost::vertex(i,graph) ] != "0";
		}
		  // report solution
		if (allFilled) {
			std::cout <<"\n Found a solution:\n\n";
			printGrid(std::cout, graph, 9,9);
			std::cout <<std::endl;
			  // only one solution wanted
			return SOLUTION_STOP;
		}
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
	InfoMode infoMode = OUT_NORMAL;
	
	std::istream* inGrid = &std::cin;
	std::ifstream* inGridFile = NULL;
	
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
		if (opts.getBoolVal("v")) {
			infoMode = OUT_VERBOSE;
		}
	} else {
		return -1;
	}
	
	try {


		//////////////////////////////////////////////////////////////
		// set up input and output streams
		//////////////////////////////////////////////////////////////

		  // set GRID input stream
		if (opts.getStrVal("sudoku").size() == 0) {
			throw ArgException("no input grid file given");
		} else if (opts.getStrVal("sudoku").compare("STDIN") != 0) {
			inGridFile = new std::ifstream(	opts.getStrVal("sudoku").c_str()
										, std::ifstream::in );
			if (!inGridFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open input grid file '" <<opts.getStrVal("sudoku") <<"'";
				throw ArgException(oss.str());
			}
			inGrid = inGridFile;
		} else {
			  // read graphs from STDIN
			inGrid = &std::cin;
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
//		sgm::Graph_Interface::CompLabel compLabel;
//		for (size_t r=0; r<rules.size(); ++r) {
//			ggl::LeftSidePattern* pattern = new ggl::LeftSidePattern(rules[r]);
//			size_t compNumber = sgm::Graph_Interface::connectedComponents( *pattern, compLabel );
//			if (rulePattern.size() <= compNumber) {
//				rulePattern.resize(compNumber+1);
//			}
//			rulePattern[compNumber].push_back(pattern);
//		}

		Graph_container c1;
		Graph_container c2;
		Graph_container& targetGraphs = c1;
		Graph_container& producedGraphs = c2;

		size_t numRow = 0, numCol = 0;

		  // parse the input grid
		GRAPHTYPE inputSudoku;
		parseGrid( *inGrid, inputSudoku, numRow, numCol );

		if (numRow != 9 || numCol != 9) {
			throw ArgException("ERROR : input Sudoku has to be a 9x9 table of [0-9]!");
		}

		  // get number of iterations == number of 0 entries
		size_t iterations = 0;
		{
			boost::property_map< GRAPHTYPE, ggl::PropNodeLabel >::type nodeLabel = boost::get( ggl::PropNodeLabel(), inputSudoku );
			for (size_t i=0; i<81; ++i ) {
				if (nodeLabel[ boost::vertex(i,inputSudoku) ] == "0") {
					iterations++;
				}
			}
		}
		  // check if nothing to do
		if (iterations == 0) {
			throw ArgException("ERROR : no unassigned fields (label == 0) within input Sudoku!");
		}

		  // store input as first target
		targetGraphs[GS_GridString::getGridString( inputSudoku )] = inputSudoku;

		*out <<"\n Input Sudoku \n\n";
		printGrid( *out, inputSudoku, 9, 9 );
		*out <<std::endl;

		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////


		  // set up graph matcher
		sgm::SGM_vf2 sgm;


		if (!opts.argExist("bfs")) {

			Sudoku_DFS::SearchTrace dfs_trace;

			Sudoku_DFS catalan_dfs;
			SudokuVisitor visitor;

			 // setup dummy graph storage for the solution
			class DoNothing : public ggl::Graph_Storage {
			public: virtual void add( const GRAPHTYPE & graph ) {};
			} solutionStorage;

			  // perform DFS
			bool solutionFound = catalan_dfs.findSolution( rules, inputSudoku, solutionStorage, visitor, &dfs_trace, true, &sgm );

			if ( !solutionFound ) {
				*out <<" no solution found ... \n"<<std::endl;
			}

		} else {

		  // print initial grid
		(*out) <<"\n";
		printGrid( *out, inputSudoku, numRow, numCol );
		(*out) <<std::endl;

		  // the graph matcher to use to unify the rule application results
		sgm::GM_vf2 graphMatcher;

		(*out) <<" " <<iterations<<" empty fields to start with \n" <<std::endl;

		  // perform iterations
		size_t it = 0;
		for (it=0; it<iterations; ++it ) {

			  // clear produced graphs container
			producedGraphs.clear();
			  // set up storage interface for MR_ApplyRule
			GS_GridString gs( producedGraphs );

			  // set up Rule applyer
			ggl::MR_ApplyRule mr_applyRule( gs );

			  // FOR EACH set of rules with equal number of connected components
			  //          in leftsidepattern :
			  // find matches and apply rules
			GraphP_container targetGraphs2(targetGraphs.size(),NULL);
			size_t i=0;
			for (Graph_container::const_iterator t=targetGraphs.begin(); t!=targetGraphs.end(); ++t,++i) {
				targetGraphs2[i] = &(t->second);
			}
			  // apply rules
			for (Graph_container::const_iterator t=targetGraphs.begin(); t!=targetGraphs.end(); ++t) {
				sgm::Graph_boost< GRAPHTYPE > target(t->second);
				  // apply rules
				sgm.findMatches( rulePattern, target, mr_applyRule, UINT_MAX );
			}

			if (infoMode >= OUT_NORMAL) {
				  // output iteration count
				(*out) <<" "<<producedGraphs.size() <<" resulting graphs with "
						<<(iterations-it-1)<<" empty fields" <<std::endl;
//				printGrid( *out, nextGrid, numRow, numCol );
				(*out)	<<std::endl;
			}

			if (producedGraphs.size() == 0) {
				if (infoMode >= OUT_NORMAL) {
					(*out) <<"\n stopping iteration since no further valid labeling possible ...\n" <<std::endl;
				}
				break;
			} else {
				*out <<"\n";
				printGrid(*out, producedGraphs.begin()->second, numRow, numCol);
				*out <<std::endl;
			}

			  // set next grid to use
			std::swap( targetGraphs, producedGraphs );
			producedGraphs.clear();

		}

		  // final output report
		if (it == iterations) {
			if (targetGraphs.size() > 0) {
				*out <<"the following solutions have been found\n\n";
				for(Graph_container::const_iterator g=targetGraphs.begin(); g!=targetGraphs.end(); ++g) {
					*out <<"\n";
					printGrid(*out, g->second, numRow, numCol );
					*out <<std::endl;
				}
			} else {
				*out <<" no solution found ... \n"<<std::endl;
			}
		}
		} // bfs

	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	inGrid = &std::cin;
	out = &std::cout;
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
		"Solves a Sudoku via depths- or breath-first-search given a start setup.\n"
		"\n"
		"The start setup has to be defined by a tabular format."
		" Use '-inputExample' for an example.\n"
		;
	
	allowedArgs.push_back(biu::COption(	
							"sudoku", false, biu::COption::STRING,
							"The initial sudoku input file name or 'STDIN' when to read from standard input (see -inputExample for a sample)",
							"STDIN"));
	allowedArgs.push_back(biu::COption(	
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(	
							"bfs", true, biu::COption::BOOL,
							"Use Breadth-First-Search instead of Depth-First-Search"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"inputExample", true, biu::COption::BOOL,
							"Displays an example for the Sudoku input encoding"));
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
 * Create the sudoku rules
 * 
 * @param toFill the container to add the rules to
 * 
 */
void
parseRules(	std::vector<ggl::Rule> & toFill )
{
	  // for each number --> add complex rule
	for (size_t i=1; i<=9; ++i ) {
		  // create rule for input i
		std::stringstream s;
		s <<"rule [\n"
			" ruleID \"" <<i <<" complex\"\n"
			" left [\n"
			"   node [ id 0 label \"0\" ]\n"
			" ]\n"
			" right [\n"
			"   node [ id 0 label \"" <<i <<"\" ]\n"
			" ]\n";
		for (size_t c=1; c<=9; ++c) {
			s <<" constrainAdj [\n"
				"   id 0\n"
				"   op "<<(c==i?"=":">") <<"\n"
				"   count 0\n"
				"   nodeLabels [\n"
				"     label \"" <<c <<"\"\n"
				"   ]\n"
				" ]\n";
		}
		s <<"]\n";

		  // parse rule
		std::pair<ggl::Rule, int>
				ret = ggl::Rule_GMLparser::parseRule( s.str() );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in complex rule specification "  <<i
				<<" at rule string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store rule
		toFill.push_back(ret.first);
	}

	  // for each number --> add set rule
	for (size_t i=1; i<=9; ++i ) {
		  // create rule for input i
		std::stringstream s;
		s <<"rule [\n"
			" ruleID \"" <<i <<" set\"\n"
			" left [\n"
			"   node [ id 0 label \"0\" ]\n"
			" ]\n"
			" right [\n"
			"   node [ id 0 label \"" <<i <<"\" ]\n"
			" ]\n"
			" constrainAdj [\n"
			"   id 0\n"
			"   op =\n"
			"   count 0\n"
			"   nodeLabels [\n"
			"     label \"" <<i <<"\"\n"
			"   ]\n"
			" ]\n"
			" constrainAdj [\n"
			"   id 0\n"
			"   op >\n"
			"   count 9\n"
			"   nodeLabels [\n";
		for (size_t c=1; c<=9; ++c) {
			if (c!=i)
				s <<"     label \"" <<c <<"\"\n";
		}
		s <<"   ]\n"
			" ]\n"
			"]\n";

		  // parse rule
		std::pair<ggl::Rule, int>
				ret = ggl::Rule_GMLparser::parseRule( s.str() );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in complex rule specification "  <<i
				<<" at rule string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store rule
		toFill.push_back(ret.first);
	}

	  // for each number --> add simple rule
	for (size_t i=1; i<=9; ++i ) {
		  // create rule for input i
		std::stringstream s;
		s <<"rule [\n"
			" ruleID \"" <<i <<"\"\n"
			" left [\n"
			"   node [ id 0 label \"0\" ]\n"
			" ]\n"
			" right [\n"
			"   node [ id 0 label \"" <<i <<"\" ]\n"
			" ]\n"
			" constrainAdj [\n"
			"   id 0\n"
			"   op =\n"
			"   count 0\n"
			"   nodeLabels [\n"
			"     label \"" <<i <<"\"\n"
			"   ]\n"
			" ]\n"
			"]\n";

		  // parse rule
		std::pair<ggl::Rule, int>
				ret = ggl::Rule_GMLparser::parseRule( s.str() );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in rule specification "  <<i
				<<" at rule string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store rule
		toFill.push_back(ret.first);
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
			if (nextLabel >= 0) {
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
				<<" : entry number ("<<row.size()
				<<") differs to preceding rows ("
				<<labels.rbegin()->size()<<")";
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
		for (size_t from=0; from<numCol; ++from) {
			for (size_t to=(from+1); to<numCol; ++to) {
				newEdge = boost::add_edge( boost::vertex( nodeShift+from, toFill )
										,  boost::vertex( nodeShift+to, toFill )
										, toFill );
				edgeLabel[ newEdge.first ] = "-";
			}
		}
	}
	  // add col edges
	for (size_t c=0; c<numCol; ++c) {
		for (size_t from=0; from<numCol; ++from) {
			for (size_t to=(from+1); to<numCol; ++to) {
				newEdge = boost::add_edge( boost::vertex( c+(from*numCol), toFill )
										,  boost::vertex( c+(to*numCol), toFill )
										, toFill );
				edgeLabel[ newEdge.first ] = "-";
			}
		}
	}
	  // add grid edges
	for (size_t c=0; c<3; ++c) {
		for (size_t r=0; r<3; ++r) {

			for (size_t lc1=0; lc1<3; ++lc1) { for (size_t lc2=0; lc2<3; ++lc2) {
				if (lc1 == lc2) continue;
				for (size_t lr1=0; lr1<3; ++lr1) { for (size_t lr2=0; lr2<3; ++lr2) {
					if (lr1 == lr2) continue;

					size_t idx_f = ((r*3+lr1)*numCol) + (c*3+lc1);
					size_t idx_t = ((r*3+lr2)*numCol) + (c*3+lc2);

					  // check if edge already existent
					if (boost::edge(boost::vertex( idx_f, toFill )
						,  boost::vertex( idx_t, toFill )
						, toFill ).second)
							continue;
					  // add edge
					newEdge = boost::add_edge( boost::vertex( idx_f, toFill )
											,  boost::vertex( idx_t, toFill )
											, toFill );
					edgeLabel[ newEdge.first ] = "-";
				}}
			}}
		}
	}

	  // add special labels to make grid traversable
	for (size_t r=0; r<numRow; ++r) {
		size_t nodeShift = r*numCol;
		  // set horizontal labels
		for (size_t c=1; c<numCol; ++c) {
			newEdge = boost::edge( boost::vertex( nodeShift+c-1, toFill )
									,  boost::vertex( nodeShift+c, toFill )
									, toFill );
			std::stringstream s; s<<(nodeShift+c-1);
			edgeLabel[ newEdge.first ] = s.str();
		}
		  // set vertical label in last column
		if (r+1<numRow) {
			nodeShift += numCol-1;
			newEdge = boost::edge( boost::vertex( nodeShift, toFill )
									,  boost::vertex( nodeShift+numCol, toFill )
									, toFill );
			std::stringstream s; s<<(nodeShift);
			edgeLabel[ newEdge.first ] = s.str();
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

	std::vector< std::string > labelInOrder = GS_GridString::getGridOrdered( g );

	  // print grid
	for (size_t r=0; r<numRow; ++r) {
		if (r%3 == 0) {
			out <<"+-------+-------+-------+\n";
		}
		for (size_t c=0; c<numCol; ++c) {
			if (c%3 == 0) {
				out <<"| ";
			}
			out <<(labelInOrder.at(r*numCol + c)=="0" ? " " : labelInOrder.at(r*numCol + c)) <<" ";
		}
		out <<"|\n";
	}
	out <<"+-------+-------+-------+\n";
}


//////////////////////////////////////////////////////////////////////////

void
giveInputExample()
{
	std::cout <<""
" 7 3 6   8 0 2   5 4 1 \n"
" 2 0 4   5 6 1   7 8 0 \n"
" 8 5 0   4 0 3   6 2 0 \n"
" 5 7 0   6 3 4   8 1 2 \n"
" 3 6 2   1 8 0   4 9 7 \n"
" 4 1 8   0 2 7   0 0 6 \n"
" 1 2 3   7 4 0   9 6 5 \n"
" 9 0 7   2 0 6   1 0 8 \n"
" 6 0 5   3 1 9   2 7 0 \n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////


