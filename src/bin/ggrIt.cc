
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
#include "sgm/Graph_boostV_p.hh"

#include "ggl/Rule.hh"
#include "ggl/RuleGraph.hh"
#include "ggl/MR_ApplyRule.hh"
#include "ggl/GS_STL.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/Graph_GMLparser.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/GS_SMILES.hh"

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

typedef std::vector< std::string > SMILES_container;
typedef std::back_insert_iterator< SMILES_container > SMILES_inserter;

typedef ggl::chem::Molecule GRAPHTYPE;

typedef std::vector< GRAPHTYPE > Graph_container;
typedef std::vector< const GRAPHTYPE * > GraphP_container;

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
parseGraphs(	std::istream & in
				, Graph_container & toFill
				, const size_t linesRead ) throw(std::exception);  

//////////////////////////////////////////////////////////////////////////

void
parseSMILES(	std::istream & in
				, SMILES_inserter & toFill
				, const size_t linesRead ) throw(std::exception);  

//////////////////////////////////////////////////////////////////////////

void
parseMolecules(	SMILES_container & smiles
				, std::vector<ggl::chem::Molecule> & toFill 
				) throw(std::exception);  

//////////////////////////////////////////////////////////////////////////

void
printRule( std::ostream& out, ggl::Rule & rule );

//////////////////////////////////////////////////////////////////////////

void
printRules( std::ostream& out, std::vector<ggl::Rule> & rules );

//////////////////////////////////////////////////////////////////////////

void
printGraphs( std::ostream& out, const Graph_container & graphs );

//////////////////////////////////////////////////////////////////////////

void
printGraphsGML( std::ostream& out
				, const Graph_container & graphs
				, const bool withSpaces );

//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample();

//////////////////////////////////////////////////////////////////////////

void
giveGraphGMLExample();

//////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv ) {
	using namespace std;
	
	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////
	int exitValue = 0;
	
	enum InfoMode {OUT_SILENT, OUT_NORMAL, OUT_VERBOSE};
	InfoMode infoMode = OUT_NORMAL;
	
	enum PrintMode {PRINT_GML_NORMAL, PRINT_GML_LINE, PRINT_ADJLIST};
	PrintMode printMode = PRINT_GML_NORMAL;
	
	std::istream* inRules = NULL;
	std::ifstream* inRulesFile = NULL;
	size_t inRulesLine = 0;
	std::istream* inGraphs = &std::cin;
	std::ifstream* inGraphsFile = NULL;
	size_t inGraphsLine = 0;
	
	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;
	
	size_t iterations = 0;
	 // rule pattern for each number of connected components
	std::vector< std::vector<sgm::Graph_Interface*> > rulePattern;
	
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
	if (opts.argExist("graphExample")) {
		giveGraphGMLExample();
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
		
		const bool allowAllIntra = opts.getBoolVal("allowAllIntra");
		
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
		  // set GRAPHS input stream
		if (opts.getStrVal("graphs").size() == 0) {
			throw ArgException("no target graph input file given");
		} else if (opts.getStrVal("graphs").compare("STDIN") != 0) {
			inGraphsFile = new std::ifstream(	opts.getStrVal("graphs").c_str()
										, std::ifstream::in );
			if (!inGraphsFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open target graph input file '" <<opts.getStrVal("graphs") <<"'";
				throw ArgException(oss.str());
			}
			inGraphs = inGraphsFile;
		} else {
			  // read graphs from STDIN
			inGraphs = &std::cin;
			  // check if rules already set to STDIN
			if (inGraphs == inRules) {
				std::ostringstream oss;
				oss	<<"cannot read both, RULES and GRAPHS, from STDIN";
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
		// setup output details
		//////////////////////////////////////////////////////////////
		
		switch(opts.getCharVal("outMode")) {
		case 'g' :
		case 'G' : printMode = PRINT_GML_NORMAL; break;
		case 's' :
		case 'S' : printMode = PRINT_GML_LINE; break;
		case 'a' :
		case 'A' : printMode = PRINT_ADJLIST; break;
		default:
				std::ostringstream oss;
				oss	<<"output mode outmode='" <<opts.getCharVal("outMode") <<"'"
					<<" not supported";
				throw ArgException(oss.str());
		}
		
		//////////////////////////////////////////////////////////////
		// parse Rules and Graphs from input
		//////////////////////////////////////////////////////////////
		
		std::vector<ggl::Rule> rules;
		parseRules( *inRules, rules, inRulesLine );
		sgm::Graph_Interface::CompLabel compLabel;
		for (size_t r=0; r<rules.size(); ++r) {
			ggl::LeftSidePattern* pattern = new ggl::LeftSidePattern(rules[r]);
			size_t compNumber = sgm::Graph_Interface::connectedComponents( *pattern, compLabel );
			if (rulePattern.size() <= compNumber) {
				rulePattern.resize(compNumber+1);
			}
			rulePattern[compNumber].push_back(pattern);
		}
		
		Graph_container c1;
		Graph_container c2;
		Graph_container& targetGraphs = c1;
		Graph_container& producedGraphs = c2;
		
		parseGraphs( *inGraphs, targetGraphs, inGraphsLine );
		
//		{
//			SMILES_container c1;
//			SMILES_container& targetSmiles = c1;
//	
//			SMILES_inserter targetInsert(targetSmiles,targetSmiles.end()); // for insert_iterator
//			//	SMILES_inserter insert(targetSmiles);	// for front_* or back_insert_iterator
//			
//			parseSMILES( *inGraphs, targetInsert, inGraphsLine );
//	
//			parseMolecules( targetSmiles, targetGraphs );
//			assert(targetSmiles.size() == targetGraphs.size());
//			
//		}
		
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n ######## PARSED RULES #######\n";
			printRules( *out, rules );
			(*out) <<"\n ######## PARSED GRAPHS ######\n\n";
			printGraphs( *out, targetGraphs );
			(*out) <<std::endl;
		}
		
		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////

		  // the graph matcher to use to unify the rule application results
		sgm::GM_vf2 graphMatcher;
		
		  // perform iterations
		for (size_t it=0; it<iterations; ++it ) {
			  // check if something to do
			if (targetGraphs.size() == 0)
				break;
			
			std::cout	<<"\n " <<(it) <<". iteration : target graphs = " 
						<<targetGraphs.size() <<std::endl;
			
			  // set up result container
			producedGraphs.clear();
			  // set up storage interface for MR_ApplyRule
			ggl::GS_STL_pushUnique gs( producedGraphs, graphMatcher);

			  // set up graph matcher
			sgm::SGM_vf2 sgm;
			
			  // FOR EACH set of rules with equal number of connected components
			  //          in leftsidepattern :
			  // find matches and apply rules
			for (size_t compNbr=1; compNbr<rulePattern.size(); ++compNbr) {
				  // check if there are rules with this number of components
				if (rulePattern[compNbr].size() == 0)
					continue;
					
				  // set up Rule applyer 
				  // check if multicomponent rule can be applied to a single molecule or not
				ggl::MR_ApplyRule mr_applyRule( gs, allowAllIntra);
//				ggl::MR_ApplyRule mr_applyRule( gs, (allowAllIntra ? 1 : compNbr));
				
				std::vector<sgm::Match_Reporter*> mr(rulePattern[compNbr].size(), &mr_applyRule);
				
				if (compNbr > 2) {
					const size_t numOfDuplicates = compNbr-1;
					  // set up graph wrapper to represent all molecules as one graph
					  // each molecule has to be present 'numOfDuplicates' times to allow
					  // matches of the rule between multiple molecules of one type
					  // and another one
					GraphP_container targetGraphs2(targetGraphs.size()*numOfDuplicates,NULL);
					for (size_t i=0; i<targetGraphs.size(); ++i ) {
						for (size_t s=0; s<numOfDuplicates; ++s) {
							targetGraphs2[(i*numOfDuplicates)+s] = &targetGraphs[i];
						}
					}
					sgm::Graph_boostV_p< GRAPHTYPE > allTargets(targetGraphs2);
					
					sgm.findMatches( rulePattern[compNbr], allTargets, mr, UINT_MAX );
				} else {
					GraphP_container targetGraphs2(targetGraphs.size(),NULL);
					for (size_t i=0; i<targetGraphs.size(); ++i ) {
						targetGraphs2[i] = &targetGraphs[i];
					}
					 // each molecule is represented only once
					sgm::Graph_boostV_p< GRAPHTYPE > allTargets(targetGraphs2);
					
					sgm.findMatches( rulePattern[compNbr], allTargets, mr, UINT_MAX );
				}
				 // check for rule application between molecules of one type only
				if (compNbr > 1) {
					  // set up graph wrapper to represent all molecules as one graph
					  // each molecule has to be present 'numOfDuplicates' times to allow
					  // matches of the rule between multiple molecules of one type
					  // and another one
					GraphP_container targetGraphs2(compNbr,NULL);
					
					for (size_t i=0; i<targetGraphs.size(); ++i ) {
						for (size_t s=0; s<targetGraphs2.size(); ++s) {
							targetGraphs2[s] = &targetGraphs[i];
						}
						sgm::Graph_boostV_p< GRAPHTYPE > allTargets(targetGraphs2);
						
						sgm.findMatches( rulePattern[compNbr], allTargets, mr, UINT_MAX );
					}
				}

			}
			
			  // make all produced SMILES target moleculed for the next iteration
			for (	Graph_container::const_iterator it=targetGraphs.begin();
					it != targetGraphs.end(); ++it )
			{
				gs.add(*it);
			}
			std::swap(producedGraphs, targetGraphs);

		}
		
		//////////////////////////////////////////////////////////////
		// write output to stream
		//////////////////////////////////////////////////////////////
		
		
		std::cout	<<"\n " <<(iterations) <<". iteration : target graphs = "
					<<targetGraphs.size() <<"\n" 
					<<std::endl;
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n ######## FINAL GRAPHS ######\n\n";
		}
		switch (printMode) {
			case PRINT_GML_NORMAL:
				printGraphsGML( *out, targetGraphs, true );
				break;
			case PRINT_GML_LINE:
				printGraphsGML( *out, targetGraphs, false );
				break;
			case PRINT_ADJLIST:
				printGraphs( *out, targetGraphs );
				break;
			default:
				break;
		}
		(*out)	<<std::endl;
	
	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}
	
	
	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	inRules = &std::cin;
	inGraphs = &std::cin;
	out = &std::cout;
	if (inRulesFile != NULL)		delete inRulesFile;
	if (inGraphsFile != NULL)		delete inGraphsFile;
	if (outFile != NULL)			delete outFile;
	for (size_t i=0; i<rulePattern.size(); ++i ) {
		for (size_t p=0; p<rulePattern[i].size(); ++p)
			delete rulePattern[i][p];
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
		"Reads a list of target graphs and graph grammar rules and "
		"iteratively applies the rules to generate new graphs.\n"
		"\n"
		"Both, rule and target input, has to be in GML format."
		" Use '-graphExample' or '-ruleExample' for an example. The rules and"
		" graphs in the input have to be separated via an ID line starting with"
		" a '#'.\n"
		;
	
	allowedArgs.push_back(biu::COption(	
							"rules", true, biu::COption::STRING, 
							"Graph grammar rules input file name or 'STDIN' when to read from standard input (use -ruleExample for a sample)"));
	allowedArgs.push_back(biu::COption(	
							"graphs", true, biu::COption::STRING, 
							"Target graphs input file name or 'STDIN' when to read from standard input (see -graphExample for a sample)",
							"STDIN"));
	allowedArgs.push_back(biu::COption(	
							"iter", true, biu::COption::INT, 
							"Number of rule application iterations",
							"1"));
	allowedArgs.push_back(biu::COption(	
							"allowAllIntra", true, biu::COption::BOOL, 
							"If present, all graph-internal rule applications are allowed, i.e. the application of rules with 2 or more unconnected components in the left side patter can applied to a single graph, otherwise NOT."));
	allowedArgs.push_back(biu::COption(	
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(	
							"outMode", true, biu::COption::CHAR, 
							"Output mode : readable (G)ML, (S)ingle line GML or (A)djacency list graph representations",
							"G"));
	allowedArgs.push_back(biu::COption(	
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"ruleExample", true, biu::COption::BOOL, 
							"Displays an example for the graph grammar rule GML encoding"));
	allowedArgs.push_back(biu::COption(	
							"graphExample", true, biu::COption::BOOL, 
							"Displays an example for the source graph GML encoding"));
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

/*!
 * Parses the given input stream for graphs in GML format and adds each graph
 * to the given container.
 * 
 * Each graph has to be headed by a comment line starting with a '#' character 
 * that specifies the name of the graph.
 * 
 * @param in the input stream to read from
 * @param toFill the container to add the graphs to
 * @param linesRead the number of lines already read from input (needed for 
 *                  error reporting)
 * 
 */
void
parseGraphs(	std::istream & in
				, Graph_container & toFill
				, const size_t linesRead = 0 )  throw(std::exception)
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
				toFill.push_back(ret.first);
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
parseMolecules(	SMILES_container & smiles
				, std::vector<ggl::chem::Molecule> & toFill 
				) throw(std::exception)
{
	for (	SMILES_container::const_iterator s = smiles.begin(); 
			s!= smiles.end(); ++s )
	{
		  // parse SMILES to graph
		std::pair<ggl::chem::Molecule,int> result 
			= ggl::chem::SMILESparser::parseSMILES(*s);
		  // check parsing result
		if (result.second != -1) {
			std::ostringstream oss;
			oss	<<"parsing error in SMILES string '"
				<<*s
				<<"' at position " <<result.second;
			throw ArgException(oss.str());
		}
		  // store parsed Molecule graph
		toFill.push_back(result.first);
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
printGraphs( std::ostream& out, const Graph_container & graphs )
{
	for (Graph_container::const_iterator it = graphs.begin(); it!=graphs.end(); ++it )
	{
		sgm::Graph_boost<GRAPHTYPE> g(*it);
		for (size_t i=0; i<g.getNodeNumber(); ++i ) {
			out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
					eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
			{
				out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			out <<" |\n";
		}
		out <<"\n";
	}
	out <<std::endl;
}

//////////////////////////////////////////////////////////////////////////

void
printGraphsGML( std::ostream& out, const Graph_container & graphs, const bool withSpaces )
{
	size_t i=0;
	for (Graph_container::const_iterator it = graphs.begin(); it!=graphs.end(); ++it )
	{
		out <<(i==0?"":"\n")<<"# result graph " <<i <<"\n";
		ggl::Graph_GML_writerT<GRAPHTYPE>::write(out, *it, withSpaces);
		i++;
	}
	out <<std::endl;
}

//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample()
{
	std::cout <<"\n\
====== RULE TO ENCODE =======================\n\
\n\
  Graph grammar rule :\n\
\n\
  0(C)   1(C) = 2(C)               0(C) - 1(C) - 2(C) \n\
                                                      \n\
   ||                     ==>       |             |   \n\
                                                      \n\
  5(C) - 4(C) = 3(C)               5(C) = 4(C) - 3(C) \n\
\n\
======= RULE IN GML =========================\n\
\n\
# My graph grammar rule ID\n\
rule [\n\
 context [\n\
   node [ id 0 label \"C\" ]\n\
   node [ id 1 label \"C\" ]\n\
   node [ id 2 label \"C\" ]\n\
   node [ id 3 label \"C\" ]\n\
   node [ id 4 label \"C\" ]\n\
   node [ id 5 label \"C\" ]\n\
 ]\n\
 left [\n\
   edge [ source 0 target 1 label \"=\" ]\n\
   edge [ source 1 target 2 label \"-\" ]\n\
   edge [ source 2 target 3 label \"=\" ]\n\
   edge [ source 4 target 5 label \"=\" ]\n\
 ]\n\
 right [\n\
   edge [ source 0 target 1 label \"-\" ]\n\
   edge [ source 1 target 2 label \"=\" ]\n\
   edge [ source 2 target 3 label \"-\" ]\n\
   edge [ source 3 target 4 label \"-\" ]\n\
   edge [ source 4 target 5 label \"-\" ]\n\
   edge [ source 5 target 0 label \"-\" ]\n\
 ]\n\
]\n\
\n\
=============================================\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////

void
giveGraphGMLExample()
{
	std::cout <<"\n\
====== GRAPH TO ENCODE =======================\n\
\n\
  Target graph :\n\
\n\
  0(A) -1- 1(B) -1- 2(A)\n\
   |        |      /    \n\
   2        2    2      \n\
   |        |  /        \n\
  3(C) -1- 4(C)         \n\
\n\
======= GRAPH IN GML =========================\n\
\n\
# My target graph ID \n\
graph [\n\
  node [ id 0 label \"A\" ]\n\
  node [ id 1 label \"B\" ]\n\
  node [ id 2 label \"A\" ]\n\
  node [ id 3 label \"C\" ]\n\
  node [ id 4 label \"C\" ]\n\
  edge [ source 0 target 1 label \"1\" ]\n\
  edge [ source 0 target 3 label \"2\" ]\n\
  edge [ source 1 target 2 label \"1\" ]\n\
  edge [ source 1 target 4 label \"2\" ]\n\
  edge [ source 2 target 4 label \"2\" ]\n\
  edge [ source 3 target 4 label \"1\" ]\n\
]\n\
\n\
==============================================\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////


