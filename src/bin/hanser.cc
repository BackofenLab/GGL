
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "biu/OptionParser.hh"
#include "biu/Timer.hh"

#include "sgm/GM_vf2.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"

#include "ggl/Rule.hh"
#include "ggl/RuleGraph.hh"
#include "ggl/MR_ApplyRule.hh"
#include "ggl/GS_STL.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/Graph_GMLparser.hh"

#include <boost/lexical_cast.hpp>

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
typedef std::vector< const sgm::Pattern_Interface * > RulePatternMap;
  //////////////////////////////////////////////////////////////////////////

typedef ggl::Graph GRAPHTYPE;

typedef std::map< std::string, GRAPHTYPE > Graph_container;
typedef std::vector< const GRAPHTYPE * > GraphP_container;

//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );  

//////////////////////////////////////////////////////////////////////////

void
parseRules(	std::vector<ggl::Rule *> & toFill );

//////////////////////////////////////////////////////////////////////////

void
parseRules(	std::istream & in, std::vector<ggl::Rule *> & toFill );

//////////////////////////////////////////////////////////////////////////

void
parseGraph(	std::istream & in
			, GRAPHTYPE & toFill
			, const size_t & linesRead = 0 ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

void
createGraphP( const GRAPHTYPE  &  graph, GRAPHTYPE & toFill );

//////////////////////////////////////////////////////////////////////////

void
createCrossingRule( const size_t crossingSize
					, const bool selfloop
					, const bool transformVertex
					, std::vector<ggl::Rule *> & rules
					, RulePatternMap & rulePattern
					, const int insertPos = -1 );

//////////////////////////////////////////////////////////////////////////

void
printRule( std::ostream& out, ggl::Rule & rule );

//////////////////////////////////////////////////////////////////////////

void
printRules( std::ostream& out, std::vector<ggl::Rule *> & rules );

//////////////////////////////////////////////////////////////////////////

void
printGraphGML(	std::ostream & out
			, const GRAPHTYPE & g
			, const bool whitespaces );

//////////////////////////////////////////////////////////////////////////

void
giveInputExample();

//////////////////////////////////////////////////////////////////////////

class GS_GraphCopy : public ggl::Graph_Storage {

	GRAPHTYPE & toFill;
	bool graphCopied;
public:
	GS_GraphCopy( GRAPHTYPE& toFill )
	 :	toFill(toFill), graphCopied(false)
	{}

	bool
	graphWasCopied(void) const {
		return graphCopied;
	}

	void
	add( const GRAPHTYPE & graph )
	{
		  // copy reported graph
		toFill = graph;
		  // note that a graph was stored
		graphCopied = true;
	}
};


//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv ) {
	using namespace std;
	
	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////
	int exitValue = 0;
	
	enum InfoMode {OUT_SILENT, OUT_NORMAL, OUT_VERBOSE};
	InfoMode infoMode = OUT_NORMAL;
	
	std::istream* inGraph = &std::cin;
	std::ifstream* inGraphFile = NULL;
	
	std::istream* inRules = &std::cin;
	std::ifstream* inRulesFile = NULL;

	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;
	
	 // rules and according patterns
	std::vector<ggl::Rule *> rules;
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
		if (opts.getStrVal("graph").size() == 0) {
			throw ArgException("no input grid file given");
		} else if (opts.getStrVal("graph").compare("STDIN") != 0) {
			inGraphFile = new std::ifstream(	opts.getStrVal("graph").c_str()
										, std::ifstream::in );
			if (!inGraphFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open input grid file '" <<opts.getStrVal("graph") <<"'";
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

		if (opts.argExist("rules")) {
			// parse rules from stream
			  // set rules input stream
			if (opts.getStrVal("rules").size() == 0) {
				throw ArgException("no input rules file given");
			} else if (opts.getStrVal("rules").compare("STDIN") != 0) {
				inRulesFile = new std::ifstream(	opts.getStrVal("rules").c_str()
											, std::ifstream::in );
				if (!inRulesFile->is_open()) {
					std::ostringstream oss;
					oss	<<"cannot open input rules file '" <<opts.getStrVal("rules") <<"'";
					throw ArgException(oss.str());
				}
				inRules = inRulesFile;
			} else {
				  // read graphs from STDIN
				inRules = &std::cin;
			}
			if (inRules == inGraph) {
				throw ArgException("cannot read graph and rules from same stream");
			}
			parseRules( *inRules, rules );
		} else {
			  // read predefined rules
			parseRules( rules );
		}

		  // generate left side pattern of each rule needed for its application
		for (size_t r=0; r<rules.size(); ++r) {
			ggl::LeftSidePattern* pattern = new ggl::LeftSidePattern(*(rules.at(r)));
			  // store in the pattern list according to the component number
			rulePattern.push_back( pattern );
		}

		  // hold the maximal crossing size for that rules are available
		size_t maxCrossingRules = 3;
		int endOfMergeRules = rules.size();
		  // create first vertex crossing transform rules
		createCrossingRule( 2, false, true, rules, rulePattern );
		createCrossingRule( 3, false, true, rules, rulePattern );
		  // create first edge crossing breaking rules
		createCrossingRule( 3, false, false, rules, rulePattern, endOfMergeRules );
		++endOfMergeRules;
		createCrossingRule( 3, true, false, rules, rulePattern, endOfMergeRules );
		++endOfMergeRules;

		  // print rules if needed
		if (infoMode == OUT_VERBOSE) {
			printRules(*out, rules);
		}

		  // parse the input graph
		GRAPHTYPE inputGraph;
		parseGraph( *inGraph, inputGraph );
		*out <<"\n Input Graph \n\n";
		printGraphGML( *out, inputGraph, true);

		  // setup graphs for processing
		GRAPHTYPE g1, g2;
		GRAPHTYPE *curGraph = &g1, *nextGraph = &g2;
		  // create initial P graph
		createGraphP( inputGraph, *curGraph );
		  // report P-graph if needed
		if( infoMode == OUT_VERBOSE) {
			*out <<"\n Initial P Graph \n\n";
			printGraphGML( *out, *curGraph, true);
		}

		*out <<std::endl;

		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////

		typedef std::set<std::string> RingLabels;
		typedef std::set< RingLabels > RingContainer;

		  // container that will hold the identified rings in the end
		RingContainer rings;

		  // set up graph matcher
		sgm::SGM_vf2 sgm;


		  // general loop that enables the creation of additional rules if needed
		while( boost::num_vertices(*curGraph) > 1 ) {
			// perform rule application
			bool operationPerformed = false;
			do {
				operationPerformed = false;

				boost::property_map<GRAPHTYPE, ggl::PropNodeLabel >::type
					curNodeLabel = boost::get( ggl::PropNodeLabel(), *curGraph );
				boost::property_map<GRAPHTYPE, ggl::PropEdgeLabel >::type
					curEdgeLabel = boost::get( ggl::PropEdgeLabel(), *curGraph );
				GRAPHTYPE::vertex_iterator curNode, nodeEnd;
				GRAPHTYPE::out_edge_iterator curEdge, edgeEnd;

				  // prune obsolete edges
				bool edgeDeleted = false;
				bool somethingPruned = false;

				// prune obsolete edges
				// i.e. prune E connected to only two adjacent vertices with the same edge label
				// i.e. prune E with non-ring loops (at least three equal node labels in the edge labels)
				do {
					edgeDeleted = false;
					for (boost::tie(curNode,nodeEnd) = boost::vertices( *curGraph ); !edgeDeleted && curNode != nodeEnd; ++curNode) {
						  // check if ring node
						if (curNodeLabel[*curNode] == "ring") {
							  // get all ring edges (duplicates due to out edge iteration)
							RingLabels ringLabel;
							for (boost::tie(curEdge,edgeEnd)=boost::out_edges( *curNode,*curGraph); curEdge!=edgeEnd; ++curEdge) {
								  // copy loop labels
								if (boost::source(*curEdge,*curGraph)==boost::target(*curEdge,*curGraph)) {
									ringLabel.insert(curEdgeLabel[*curEdge]);
								}
							}
							  // add ring to container
							rings.insert(ringLabel);
							  // delete ring node
							while( boost::edge(*curNode,*curNode,*curGraph).second ) {
								boost::remove_edge(*curNode,*curNode,*curGraph);
							}
							boost::clear_vertex( *curNode, *curGraph );
							boost::remove_vertex( *curNode, *curGraph );
							edgeDeleted = true;
							somethingPruned = true;
							break;
						} else
							// check if an edge node to prune
						if (curNodeLabel[*curNode] == "E" ) {
							  // get out edge information to adjacent nodes (ignore self loops)
							size_t realOutEdges = 0;
							std::set< std::string > realOutLabel;
							std::map< std::string, double > nodeCount;
							for(boost::tie(curEdge,edgeEnd)=boost::out_edges( *curNode,*curGraph); curEdge!=edgeEnd; ++curEdge) {
								const std::string & label = curEdgeLabel[*curEdge];
								if (boost::target(*curEdge,*curGraph) != boost::source(*curEdge,*curGraph)) {
									realOutEdges += 1;
									realOutLabel.insert(label);
								}
								size_t cutPos = label.find('-');
								  // store first node index
								std::string nodeLabel = label.substr(0,cutPos);
								  // check if node already known, otherwise add to container
								if (nodeCount.find(nodeLabel)==nodeCount.end()) {
									nodeCount[nodeLabel] = 0;
								}
								  // increase counter
								if (boost::target(*curEdge,*curGraph) != boost::source(*curEdge,*curGraph)) {
									nodeCount[nodeLabel] += 1;
								} else {
									nodeCount[nodeLabel] += 0.5;
								}
								  // store second node index
								nodeLabel = label.substr(cutPos+1,label.size());
								  // check if node already known, otherwise add to container
								if (nodeCount.find(nodeLabel)==nodeCount.end()) {
									nodeCount[nodeLabel] = 0;
								}
								  // increase counter
								if (boost::target(*curEdge,*curGraph) != boost::source(*curEdge,*curGraph)) {
									nodeCount[nodeLabel] += 1;
								} else {
									  // self loop edges will be seen twice -> count only half
									nodeCount[nodeLabel] += 0.5;
								}
							}
							  // check if path with equal edges to two different vertices
							bool deleteEdgeNode = (realOutLabel.size() == 1) && (realOutEdges == 2);
							  // check if multiple node occurrences
							if (!deleteEdgeNode && realOutEdges == 2) {
								for (std::map< std::string, double >::const_iterator nc=nodeCount.begin(); !deleteEdgeNode && nc!=nodeCount.end(); ++nc) {
									  // check if this node was seen more than twice
									deleteEdgeNode = nc->second > 2;
								}
							}
							  // check if vertex is to be deleted
							if ( deleteEdgeNode ) {
								  // clear self loops (not done by clear_vertex)
								while( boost::edge(*curNode,*curNode,*curGraph).second ) {
									boost::remove_edge(*curNode,*curNode,*curGraph);
								}
								  // clear remaining in and out edges
								boost::clear_vertex( *curNode, *curGraph );
								  // delete whole node
								boost::remove_vertex( *curNode, *curGraph );
								edgeDeleted = true;
								somethingPruned = true;
								break;
							}

						}
					}
				} while (edgeDeleted);

				if (infoMode == OUT_VERBOSE && somethingPruned) {
					  // print pruned graph
					*out <<"\n edge pruning\n"<<std::endl;
					printGraphGML( *out, *curGraph, true );
				}

				  // setup matching target
				sgm::Graph_boost<GRAPHTYPE> curTarget(*curGraph);


				/*
				 *
				 * TODO change VF2 search
				 * - MC_EdgeLabel constraints into EdgeLabelData
				 *   --> check directly within compatible(..)
				 * - add nodes in decreasing node degree / constraint number order to search graph
				 *
				 *
				 * TODO aromatic edge support in ITS
				 *
				 */

				 // check each rule until a match was found
				for (size_t i=0; !operationPerformed && i<rulePattern.size(); ++i) {

					  // setup graph storage
					GS_GraphCopy gs_copy(*nextGraph);
					  // setup rule application system
					ggl::MR_ApplyRule mr_applyRule(gs_copy);

					  // perform matching -> only one hit needed to proceed
					sgm.findMatches( *(rulePattern.at(i)), curTarget, mr_applyRule, 1);
					  // check if an operation was performed
					operationPerformed = gs_copy.graphWasCopied();
					  // report which operation was done
					if (operationPerformed && infoMode == OUT_VERBOSE) {
						*out <<"\n rule : " <<rules.at(i)->getID() <<std::endl;
						  // report resulting graph
						printGraphGML( *out, *nextGraph, true );
					}
				}
				  // make result graph the input for the next iteration
				if (operationPerformed) {
					  // swap current and next for the following iteration
					std::swap(curGraph,nextGraph);
				}

			} while (operationPerformed);

			  // check if a solution was found
			  // --> otherwise: the current rule set is not sufficient
			if ( boost::num_vertices(*curGraph) > 1 ) {

				  // add additional rules that might enable the complete solution
				maxCrossingRules +=1;
				  // create vertex crossing transform rule
				createCrossingRule( maxCrossingRules, false, true, rules, rulePattern );
				  // create edge crossing breaking rules
				createCrossingRule( maxCrossingRules, false, false, rules, rulePattern, endOfMergeRules );
				++endOfMergeRules;
				createCrossingRule( maxCrossingRules, true, false, rules, rulePattern, endOfMergeRules );
				++endOfMergeRules;
			}

		}

		  // print all rings
		*out <<"\n number of rings = " <<rings.size() <<"\n" <<std::endl;
		for (RingContainer::const_iterator r=rings.begin(); r!=rings.end(); ++r) {
			  // print all ring edge label
			*out <<" ring =";
			for (std::set<std::string>::const_iterator l=r->begin(); l!=r->end(); ++l) {
				*out <<" " <<*l;
			}
			*out <<std::endl;

		}

	} catch (std::exception& ex) {
		std::cerr <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}


	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	inRules = &std::cin;
	inGraph = &std::cin;
	out = &std::cout;
	if (inRulesFile != NULL)		delete inRulesFile;
	if (inGraphFile != NULL)		delete inGraphFile;
	if (outFile != NULL)			delete outFile;
	for (size_t p=0; p<rulePattern.size(); ++p) {
		  // first delete pattern
		if (rulePattern[p] != NULL) delete rulePattern[p];
		  // now delete rule
		if (rules[p] != NULL) delete rules[p];
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
		"Performs a ring perception based on Graph Grammar application to enumerate all rings"
		" within a given graph using a Graph Grammar Rule encoding of the algorithm by Hanser et al. (1996).\n"
		"\n"
		"The graph has to be given in GML format."
		" Use '-inputExample' for an example.\n"
		;
	
	allowedArgs.push_back(biu::COption(	
							"graph", false, biu::COption::STRING,
							"The initial graph input file name or 'STDIN' when to read from standard input (see -inputExample for a sample)",
							"STDIN"));
	allowedArgs.push_back(biu::COption(	
							"rules", true, biu::COption::STRING,
							"The rules to be used instead of the predefined ones. A file name or 'STDIN' is expected."));
	allowedArgs.push_back(biu::COption(
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"inputExample", true, biu::COption::BOOL,
							"Displays an example for the graph GML input encoding"));
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
 * helper function to generate rule strings for edge compression rules
 *
 * @param thisSelfLoop whether or not the edge node to delete should have a
 *                     self-loop
 * @param neighSelfLoop whether or not the neighbored edge node should have a
 *                     self-loop
 * @param isLoop whether or not the edge node and the neighbored edge node form
 *                     together with another "V" node a triangle
 * @param ruleStrings the container to add the rule string to
 *
 */
void
createPathCompressString( const bool thisSelfLoop
						, const bool neighSelfLoop
						, const bool isLoop
						, std::vector<std::string> & ruleStrings )
{
	std::stringstream rule;

	  // create rule string
	rule <<"rule [\n"
			<<" ruleID \"merging "
			<<(isLoop?"within triangle ring ":"")
			<<(thisSelfLoop?"self-looped ":"")
			<<"edge node with "
			<<(neighSelfLoop?"self-looped ":"")
			<<"neighbored edge node\"\n"
			<<" wildcard \"*\"\n"
			<<" context [\n"
			<<"   node [ id 1 label \""<<(isLoop?"V":"*")<<"\" ]\n"
			<<"   node [ id 2 label \"E\" ]\n"
			<<(isLoop?"":"   node [ id 3 label \"*\" ]\n")
			<<"   edge [ source 2 target "<<(isLoop?1:3)<<" label \"*\" ]\n"
			<<(neighSelfLoop?"   edge [ source 2 target 2 label \"*\" ]\n":"")
			<<" ]\n"
			<<" left [\n"
			<<"   node [ id 0 label \"E\" ]\n"
			<<(thisSelfLoop?"   edge [ source 0 target 0 label \"*\" ]\n":"")
			<<"   edge [ source 0 target 1 label \"*\" ]\n"
			<<"   edge [ source 0 target 2 label \"*\" ]\n"
			<<" ]\n"
			<<" constrainAdj [ id 0 op = count " <<(thisSelfLoop ? 3 : 2)
			<<" nodeLabels [ label \"*\" ] ]\n"
			<<" constrainAdj [ id 2 op = count " <<(neighSelfLoop ? 3 : 2)
			<<" nodeLabels [ label \"*\" ] ]\n"
			<<" copyAndPaste [ source 0 id 2 ]\n"
			<<"]\n"
			;

	  // add to rule string container
	ruleStrings.push_back(rule.str());
}





/*!
 * Create the basic Hanser ring perception rules
 * 
 * @param toFill the container to add the rules to
 * 
 */
void
parseRules(	std::vector<ggl::Rule *> & toFill )
{
	  // compress
	std::vector < std::string > ruleStrings;

	ruleStrings.push_back( std::string(
		"rule [\n"
		" ruleID \"remove single node\"\n"
		" wildcard \"*\"\n"
		" left [\n"
		"   node [ id 0 label \"V\" ]\n"
		" ]\n"
		" constrainAdj [ id 0 op = count 0 edgeLabels [ label \"*\" ] ]\n"
		"]\n"
	));

	ruleStrings.push_back( std::string(
		"rule [\n"
		" ruleID \"remove dangling end\"\n"
		" wildcard \"*\"\n"
		" left [\n"
		"   node [ id 1 label \"*\" ]\n"
		" ]\n"
		" constrainAdj [ id 1 op = count 1 edgeLabels [ label \"*\" ] ]\n"
		"]\n"
	));


	ruleStrings.push_back( std::string(
		"rule [\n"
		" ruleID \"close and decouple ring without self-loop\"\n"
		" wildcard \"*\"\n"
		" context [\n"
		"   node [ id 0 label \"V\" ]\n"
		" ]\n"
		" left [\n"
		"   node [ id 1 label \"E\" ]\n"
		"   edge [ source 0 target 1 label \"*\" ]\n"
		"   edge [ source 0 target 1 label \"*\" ]\n"
		" ]\n"
		" constrainAdj [ id 1 op = count 1 nodeLabels [ label \"*\" ] ]\n"
		" right [\n"
		"   node [ id 1 label \"ring\" ]\n"
		" ]\n"
		" copyAndPaste [ source 0 target 1 id 1 ]\n"
		"]\n"
	));

	ruleStrings.push_back( std::string(
		"rule [\n"
		" ruleID \"close and decouple ring with self-loops\"\n"
		" wildcard \"*\"\n"
		" context [\n"
		"   node [ id 0 label \"V\" ]\n"
		" ]\n"
		" left [\n"
		"   node [ id 1 label \"E\" ]\n"
		"   edge [ source 0 target 1 label \"*\" ]\n"
		"   edge [ source 0 target 1 label \"*\" ]\n"
		"   edge [ source 1 target 1 label \"*\" ]\n"
		" ]\n"
		" constrainAdj [ id 1 op = count 2 nodeLabels [ label \"*\" ] ]\n"
		" right [\n"
		"   node [ id 2 label \"ring\" ]\n"
		" ]\n"
		" copyAndPaste [ source 1 id 2 ]\n"
		" copyAndPaste [ source 0 target 1 id 2 ]\n"
		"]\n"
	));


	  // create all combinations of ring path compressions
	createPathCompressString( false, false,  true, ruleStrings );
	createPathCompressString(  true, false,  true, ruleStrings );
	createPathCompressString(  true,  true,  true, ruleStrings );

	  // create all combinations of path compressions
	createPathCompressString( false, false, false, ruleStrings );
	createPathCompressString(  true, false, false, ruleStrings );
	createPathCompressString( false,  true, false, ruleStrings );
	createPathCompressString(  true,  true, false, ruleStrings );


	  // for each rule string --> add according rule
	for (size_t i=0; i<ruleStrings.size(); ++i) {

		try {
			  // parse rule
			std::pair<ggl::Rule, int>
				ret = ggl::Rule_GMLparser::parseRule( ruleStrings.at(i) );
			  // check for parsing error
			if (ret.second != -1) {
				std::ostringstream oss;
				oss	<<" Rule parsing error in rule specification "  <<i
					<<" at rule string position " <<ret.second;
				throw ArgException(oss.str());
			}
			  // store rule
			toFill.push_back(new ggl::Rule(ret.first));
		} catch (ggl::Rule_GML_error & error) {
			std::ostringstream oss;
			oss	<<" Rule parsing error in rule specification "  <<i
				<<" : " <<error.what() <<"\n"<<ruleStrings.at(i);
			throw ArgException(oss.str());
		}
	}

}


//////////////////////////////////////////////////////////////////////////

void
parseRules(	std::istream& in, std::vector<ggl::Rule *> & toFill )
{
	std::string ruleString = std::string(""), line;
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );

		  // ignore comment lines
		if (line.size()>0 && line.at(0)=='#') {
			continue;
		}

		  // check if next rule header or end of file was reached
		if (line.find("rule [") != std::string::npos || in.eof()) {
			size_t cutPos  = line.find("rule [");
			std::string nextRuleString = "";
			if (cutPos != std::string::npos) {
				  // append remaining rule string
				ruleString.append( line.substr(0,cutPos) );
				  // store beginning of next rule string
				nextRuleString = line.substr(cutPos, line.size());
			}
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
				toFill.push_back(new ggl::Rule(ret.first));
			}
			  // clear processed rule
			ruleString.clear();
			  // set up with beginning of next rule string if present
			ruleString = nextRuleString;
		} else {
			  // append line content to rule string
			ruleString.append(line);
		}
	}

}

//////////////////////////////////////////////////////////////////////////

void
parseGraph(	std::istream & in
			, GRAPHTYPE & toFill
			, const size_t & linesRead ) throw(std::exception)
{
	std::string line;
	size_t lineNumber = linesRead;
	std::string graphString = std::string("");
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;

		  // ignore comment lines
		if (line.size()>0 && line.at(0)=='#') {
			continue;
		}

		  // append line content to rule string
		graphString.append(line);
	}
	  // convert rule string to rule object
	if (graphString.size() > 0) {
		std::pair<GRAPHTYPE, int>
				ret = ggl::Graph_GMLparser::parseGraph( graphString );
		  // check for parsing error
		if (ret.second != -1) {
			std::ostringstream oss;
			oss	<<" Graph parsing error in graph specification BEFORE line "
				<<line
				<<" at graph string position " <<ret.second;
			throw ArgException(oss.str());
		}
		  // store
		toFill = ret.first;
	}
}


//////////////////////////////////////////////////////////////////////////

void
createGraphP( const GRAPHTYPE  &  graph, GRAPHTYPE & toFill )
{
	boost::property_map<GRAPHTYPE, ggl::PropNodeIndex >::const_type
		inNodeIndex(boost::get( ggl::PropNodeIndex(), graph ));

	  //! node label access
	boost::property_map<GRAPHTYPE, ggl::PropNodeLabel >::const_type
		inNodeLabel(boost::get( ggl::PropNodeLabel(), graph ));
	boost::property_map<GRAPHTYPE, ggl::PropNodeLabel >::type
		outNodeLabel(boost::get( ggl::PropNodeLabel(), toFill ));

	boost::property_map<GRAPHTYPE, ggl::PropEdgeLabel>::type
		outEdgeLabel(boost::get( ggl::PropEdgeLabel(), toFill ));

	  // create nodes in P graph
	for (size_t i=0; i<boost::num_vertices(graph); ++i) {
		GRAPHTYPE::vertex_descriptor v = boost::add_vertex( toFill );
		outNodeLabel[v] = "V";
	}
	  // create edges in P graph
	GRAPHTYPE::edge_iterator curEdge, edgeEnd;
	for (boost::tie(curEdge,edgeEnd) = boost::edges( graph ); curEdge != edgeEnd; ++curEdge) {
		  // get edge adjancency information
		size_t from = inNodeIndex[boost::source(*curEdge,graph)]
		       , to = inNodeIndex[boost::target(*curEdge,graph)];
		  // create new edge
		std::pair< GRAPHTYPE::edge_descriptor, bool > e
			= boost::add_edge(
					  boost::vertex( from, toFill )
					, boost::vertex( to, toFill )
					, toFill );
		  // set edge label based on adjacent nodes
		outEdgeLabel[ e.first ] = boost::lexical_cast<std::string>(from) + "-" + boost::lexical_cast<std::string>(to);
	}

}

//////////////////////////////////////////////////////////////////////////

void
createCrossingRule( const size_t crossingSize
					, const bool selfloop
					, const bool transformVertex
					, std::vector<ggl::Rule *> & rules
					, RulePatternMap & rulePattern
					, const int insertPos )
{
	  // ensure minimal crossing sizes
	if (crossingSize<2)
		return;
	if (!transformVertex && crossingSize<3)
		return;

	std::stringstream ruleString("");

	ruleString <<"rule [\n";
	if (transformVertex) {
		ruleString	<<" ruleID \"transform vertex crossing " <<crossingSize
					<<" to edges\"\n";
	} else {
		ruleString	<<" ruleID \"break edge crossing " <<crossingSize
					<<(selfloop?" with":" without") <<" self-loops\"\n";
	}
	ruleString <<" wildcard \"*\"\n";

	  // create context
	std::stringstream context;
	context <<" context [\n";
	for (size_t i=1; i<=crossingSize; ++i) {
		context <<"  node [ id " <<i <<" label \"*\" ]\n";
	}
	context <<" ]\n";
	  // append to rule
	ruleString <<context.str();

	  // create left side
	std::stringstream left;
	left <<" left [\n";
	if (transformVertex) {
		left <<"  node [ id 0 label \"V\" ]\n";
	} else {
		left <<"  node [ id 0 label \"E\" ]\n";
	}
	if (selfloop) {
		left <<"  edge [ source 0 target 0 label \"*\" ]\n";
	}
	for (size_t i=1; i<=crossingSize; ++i) {
		left <<"  edge [ source 0 target " <<i <<" label \"*\" ]\n";
	}
	left <<" ]\n";
	  // append to rule
	ruleString <<left.str();

	  // create adjacency constraint
	ruleString	<<" constrainAdj [ id 0 op = count "
				<<(selfloop ? crossingSize+1 : crossingSize)
				<<" nodeLabels [ label \"*\" ] ]\n";

	  // create right side
	size_t nextEdgeID = crossingSize+1;
	std::stringstream right;
	right <<" right [\n";
	std::stringstream copyAndPaste;
	for (size_t i=1; i<crossingSize; ++i) {
		for (size_t j=i+1; j<=crossingSize; ++j) {
			  // create edge node
			right <<"  node [ id " <<nextEdgeID <<" label \"E\" ]\n";
			  // connect edge node via copyAndPaste
			if (selfloop) {
				copyAndPaste <<" copyAndPaste [ source 0 target 0 id "<<nextEdgeID<<" ]\n";
			}
			copyAndPaste <<" copyAndPaste [ source 0 target "<<i<<" id "<<nextEdgeID<<" ]\n";
			copyAndPaste <<" copyAndPaste [ source 0 target "<<j<<" id "<<nextEdgeID<<" ]\n";
			  // increase edge counter
			nextEdgeID++;
		}
	}
	right <<" ]\n";
	  // append to rule
	ruleString	<<right.str()
				<<copyAndPaste.str()
				<<"]\n";

	  // parse rule string
	std::pair<ggl::Rule, int>
			ret = ggl::Rule_GMLparser::parseRule( ruleString.str() );
	  // check for parsing error
	if (ret.second != -1) {
		std::ostringstream oss;
		oss	<<" Compress rule parsing error at rule string position " <<ret.second;
		throw ArgException(oss.str());
	}

	  // store new rule
	if (insertPos < 0)
		rules.push_back( new ggl::Rule(ret.first) );
	else {
		assert(insertPos < (int)rulePattern.size());
		rules.insert( rules.begin()+insertPos, new ggl::Rule(ret.first) );
	}

	  // create and store pattern for matching
	if (insertPos < 0)
		rulePattern.push_back( new ggl::LeftSidePattern(**(rules.rbegin())) );
	else {
		assert(insertPos <(int) rulePattern.size());
		rulePattern.insert( rulePattern.begin()+insertPos, new ggl::LeftSidePattern(*(rules.at(insertPos))) );
	}

}


//////////////////////////////////////////////////////////////////////////


void
printRules( std::ostream& out, std::vector<ggl::Rule *> & rules )
{
	for (size_t i=0; i<rules.size(); ++i ) {
		printRule(out, *(rules.at(i)));
		out <<"\n";
	}

}

void printCopyAndPaste( std::ostream & out, const ggl::Rule rule )
{
	out <<"\n number of nodes with copy-and-paste operations = " <<rule.getCopyAndPasteOperations().size() <<"\n"<<std::endl;
	for ( ggl::Rule::CopyAndPasteOperations::const_iterator n=rule.getCopyAndPasteOperations().begin(); n!=rule.getCopyAndPasteOperations().end(); ++n) {
		for (ggl::Rule::CopyAndPasteOperations::mapped_type::const_iterator cnp=n->second.begin(); cnp!=n->second.end(); ++cnp) {
			out <<" copy-and-paste : source("<<cnp->source <<") pasteID("<<cnp->pasteID<<") target(";
			if (cnp->target != (size_t)INT_MAX)
				out <<cnp->target;
			out <<") edgeLabels (";
			for (ggl::Rule::RuleCnP::LabelSet::const_iterator l=cnp->edgeLabels.begin(); l!=cnp->edgeLabels.end(); ++l)
				out <<" '"<<*l<<"'";
			out <<")"<<std::endl;
		}
	}
}


//////////////////////////////////////////////////////////////////////////
void printConstraints( std::ostream & out, const ggl::Rule & rule )
{
	out <<"\n number of constraints = " <<rule.getLeftSide().constraints.size() <<std::endl;

	typedef std::vector< sgm::Pattern_Interface::Match_Constraint* > CV;
	const CV & cv = rule.getLeftSide().constraints;
	for (CV::const_iterator c = cv.begin(); c != cv.end(); ++c ) {
		const sgm::MC_NodeLabel * cur = dynamic_cast< const sgm::MC_NodeLabel * >(*c);
		{
		if (cur != NULL) {
			out <<" + MC_NodeLabel "<<std::endl;
			out <<"   id " <<cur->constrainedNodeID <<std::endl;
			out <<"   op " <<(cur->compareType==sgm::MC_NodeLabel::ALLOWED?"=":"!") <<std::endl;
			out <<"   nodeLabels = " ;
			for (sgm::MC_NodeLabel::LabelSet::const_iterator l=cur->nodeLabels.begin(); l!=cur->nodeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;
			continue;
		}}{
		const sgm::MC_EdgeLabel * cur = dynamic_cast< const sgm::MC_EdgeLabel * >(*c);
		if (cur != NULL) {
			out <<" + MC_EdgeLabel "<<std::endl;
			out <<"   from " <<cur->constrainedFromID <<std::endl;
			out <<"   to   " <<cur->constrainedToID <<std::endl;
			out <<"   op " <<(cur->compareType==sgm::MC_EdgeLabel::ALLOWED?"=":"!") <<std::endl;
			out <<"   edgeLabels = " ;
			for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=cur->edgeLabels.begin(); l!=cur->edgeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;
			continue;
		}}{
		const sgm::MC_NoEdge * cur = dynamic_cast< const sgm::MC_NoEdge * >(*c);
		if (cur != NULL) {
			out <<" + MC_NoEdge "<<std::endl;
			out <<"   from " <<cur->constrainedFromID <<std::endl;
			out <<"   to   " <<cur->constrainedToID <<std::endl;
			continue;
		}}{
		const sgm::MC_NodeAdjacency * cur = dynamic_cast< const sgm::MC_NodeAdjacency * >(*c);
		if (cur != NULL) {
			out <<" + MC_NodeAdjacency "<<std::endl;
			out <<"   id " <<cur->constrainedNodeID <<std::endl;
			out <<"   op " ;
			switch(cur->op) {
			case sgm::MC_NodeAdjacency::MC_EQ : out <<"="; break;
			case sgm::MC_NodeAdjacency::MC_G : out <<">"; break;
			case sgm::MC_NodeAdjacency::MC_GQ : out <<">="; break;
			case sgm::MC_NodeAdjacency::MC_L : out <<"<"; break;
			case sgm::MC_NodeAdjacency::MC_LQ : out <<"<="; break;
			default : out <<"unknown operator type";
			}
			out <<std::endl;
			out <<"   count " <<cur->count <<std::endl;
			out <<"   nodeLabels = " ;
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=cur->nodeLabels.begin(); l!=cur->nodeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;
			out <<"   edgeLabels = " ;
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=cur->edgeLabels.begin(); l!=cur->edgeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;

			continue;
		}}
		out <<" !!! unknown constraint type !!! " <<std::endl;
	}

}


void
printRule( std::ostream& out, ggl::Rule & rule )
{
	out <<"\n RULE " <<rule.getID() <<std::endl;
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
	printConstraints( out, rule );
	printCopyAndPaste( out, rule );
	out <<"\n";
}


//////////////////////////////////////////////////////////////////////////

void
printGraphGML( std::ostream& out, const GRAPHTYPE & graph, const bool withSpaces )
{
	ggl::Graph_GML_writer::write(out, graph, withSpaces);
	out <<std::endl;
}

//////////////////////////////////////////////////////////////////////////

void
giveInputExample()
{
	std::cout <<""
" graph [\n"
"   node [ id 0 label \"A\" ] \n"
"   node [ id 1 label \"B\" ] \n"
"   node [ id 2 label \"C\" ] \n"
"   edge [ source 0 target 1 label \"-\" ] \n"
"   edge [ source 0 target 2 label \"-\" ] \n"
"   edge [ source 1 target 2 label \"-\" ] \n"
" ]\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////


