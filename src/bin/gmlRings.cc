

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
#include <sgm/RP_Hanser96.hh>

#include <ggl/Graph.hh>
#include <ggl/Graph_GML_writer.hh>
#include <ggl/Graph_GMLparser.hh>


#include "version.hh"

	using namespace ggl;


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


	typedef std::vector< Graph > Graph_container;

/////////////////////////////////////////////////////////////////////////////


void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );

//////////////////////////////////////////////////////////////////////////

void
parseGraphs(	std::istream & in
				, Graph_container & toFill
				, const size_t linesRead ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

/*!
 *
 * Dummy aromaticity perception to access rings to predict
 *
 */
class RR_list
	: public sgm::RingReporter
{
public:

	typedef std::vector< RingList > RingListVec;

protected:

	  //! container that will hold all rings of the current graph
	RingListVec allRings;


	  /*!
	   * Comparison of ring lists based on the ring size.
	   */
	class RingSizeLess {
	public:

		  /*!
		   * less comparison of two ring descriptors.
		   * @param x the ring to be checked if smaller
		   * @param y the ring to compare to
		   * @return true if x.size < y.size;
		   * 		false otherwise
		   */
		bool operator () ( const RingList & x
						  , const RingList & y )
		{
			return x.size() < y.size();
		}
	};


public:

	RR_list()
	 : allRings()
	{}

	virtual ~RR_list()
	{}

	void
	reportRing(	const sgm::Graph_Interface& graph
				, const RingList & ringList )
	{
		  // store description
		allRings.push_back( ringList );
	}


	RingListVec &
	getAllRings()
	{
		return allRings;
	}

	void
	pruneFusedRings( RR_list::RingListVec & rings );

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



		  // set graph input stream
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
		// parse graphs from input
		//////////////////////////////////////////////////////////////

		  // get master graph to match onto
		Graph_container parsedGraphs;
		parseGraphs( *in, parsedGraphs, inLine );
		if (parsedGraphs.size() != 1)  {
			throw ArgException("could not parse only one graph from 'graph'");
		}
		Graph & graph = *(parsedGraphs.begin());


		//////////////////////   RING ENUMERATION  /////////////////////////////////////


		  // store all rings up to the maximally predictable size
		sgm::RP_Hanser96 ringFinder;
		  // setup graph interface
		Graph_boost graphWrapper(graph);
		  // setup ring container
		RR_list reporter;
		  // run ring perception
		ringFinder.findRings( graphWrapper, reporter, maxRingSize );

		  // prune rings if needed
		if (opts.argExist("pruneFusedRings")) {
			  // prune rings to remove fused representation
			reporter.pruneFusedRings( reporter.getAllRings() );
		}

		//////////////////////   OUTPUT  /////////////////////////////////////

		const RR_list::RingListVec & allRings = reporter.getAllRings();

		for (size_t i=0; i< allRings.size(); ++i) {

			const sgm::RingReporter::RingList & nodes = allRings.at(i);

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
		"Enumerates all rings for a given graph.\n"
		;

	allowedArgs.push_back(biu::COption(
							"graph", true, biu::COption::STRING,
							"The name of the file that holds the graph encoding in GML or 'STDIN' when to read GML from standard input",
							"STDIN"));
	allowedArgs.push_back(biu::COption(
							"maxRingSize", true, biu::COption::INT,
							"The maximal ring size to be considered",
							"15"));
	allowedArgs.push_back(biu::COption(
							"pruneFusedRings", true, biu::COption::BOOL,
							"prunes the outer rings that emerge if two rings share a single edge"));
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

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////



	void
	RR_list::
	pruneFusedRings( RR_list::RingListVec & rings )
	{

		  // sort rings by size
		RingSizeLess ringSizeLess;

		std::sort( rings.begin(), rings.end(), ringSizeLess );

		  //! defines an edge; NOTE: use as ordered pairs, i.e. first <= second
		typedef std::pair< size_t, size_t > Edge;
		  //! container to maintain aromatic edges without duplicates
		typedef std::set< Edge > EdgeSet;

		typedef std::vector< EdgeSet > EdgeSetVec;

		EdgeSetVec edges(rings.size());

		  // generate edge list description for each ring
		for (size_t i=0; i<rings.size(); ++i) {
			RingList::const_iterator cur=rings.at(i).begin(), last=rings.at(i).begin();
			for (cur++; cur!=rings.at(i).end(); ++cur,++last) {
				if (*cur<*last) {
					edges[i].insert( Edge(*cur,*last) );
				} else {
					edges[i].insert( Edge(*last,*cur) );
				}
			}
		}

//std::cerr <<"\nNEXT\n";
		  // prune rings
		size_t numOfprunedRings = 0;
		for (int i=rings.size()-1; i>=0; i--) {
//std::cerr <<"next ring = ";
//for (RingList::const_iterator x=rings.at(i).begin(); x!=rings.at(i).end(); ++x) {
//	std::cerr <<*x <<" ";
//}
//std::cerr <<"\n";
			bool canBePruned = false;
			for( size_t c=0; !canBePruned && c < (size_t)i; ++c ) {
//std::cerr <<"  compare to ring = ";
//for (RingList::const_iterator x=rings.at(c).begin(); x!=rings.at(c).end(); ++x) {
//	std::cerr <<*x <<" ";
//}
//std::cerr <<"\n";
				  // skip already deleted rings
				if (edges.at(c).empty()) {
					continue;
				}
				  // stop if ring is of equal size or larger
				if (edges.at(c).size() >= edges.at(i).size()) {
					break;
				}
				  // check for containment of this ring c within ring i

				  // get iterators on edge sets to compare
				EdgeSet::const_iterator largeEdge = edges.at(i).begin(), largeEnd = edges.at(i).end();
				EdgeSet::const_iterator curEdge = edges.at(c).begin(), curEnd = edges.at(c).end();
				  // sequential check
				size_t curOverlap = 0;
				  // check only until up to two differences have been found
				while(largeEdge != largeEnd && curEdge != curEnd) {
					  // check if overlap
					if (*largeEdge == *curEdge) {
						++curOverlap;
						++largeEdge;
						++curEdge;
					  // check which counter to increase
					} else if (*largeEdge < *curEdge) {
						++largeEdge;
					} else {
						++curEdge;
					}
				}
				  // check if smaller current ring is contained excluding one edge
				  // if only one edge is left : ring is contained
				canBePruned = (curOverlap+1 == edges.at(c).size());

			}
			  // prune ring if obsolete
			if (canBePruned) {
				  // remove edge information for ring at position i
				edges[i].clear();
				  // count pruning
				numOfprunedRings++;
//std::cerr <<"==> pruned\n";
			}
		}


		  // shift remaining rings to the front
		if (numOfprunedRings>0) {
			size_t fillPos = 0;
			for( size_t i=0; i<rings.size(); ++i ) {
				if (!edges.at(i).empty()) {
					  // check if overwrite is needed
					if (fillPos < i) {
						rings[fillPos] = rings[i];
					}
					  // increase overwrite counter
					fillPos++;
				}
			}
			  // drop deleted elements
			rings.resize(fillPos);
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
