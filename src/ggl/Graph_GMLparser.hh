#ifndef GGL_GRAPH_GML_PARSER_HH_
#define GGL_GRAPH_GML_PARSER_HH_

#include <utility>
#include <vector>

#include "ggl/Graph.hh"

namespace ggl {


	/*! @brief Graph GML parser
	 *
	 * This is a wrapper for the ggl::Graph_GML_grammar BNF parser that parses
	 * a GML string representation of a ggl::Graph object. See there for
	 * further details.
	 * 
	 * @author Martin Mann (c) 2012 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */ 
	class Graph_GMLparser
	{
		
	public:
	    
	      //! Parses a GML string and generates a Graph object
	      //! @param GML_string the string to parse
	      //! @return pair.first = the graph object
	      //!         pair.second = -1 if parsing was successful,
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
	    static
	    std::pair< Graph, int >
	    parseGraph( const std::string & GML_string );
	    
	};


} // namespace ggl


#endif /*GGL_GRAPH_GML_PARSER_HH_*/
