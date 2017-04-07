
#include "ggl/Graph_GMLparser.hh"
#include "ggl/Graph_GML_grammar.hh"

namespace ggl {

//##############################################################################
		
	std::pair< Graph, int >
	Graph_GMLparser
    ::parseGraph( const std::string & GML_string )
    {
		// forward call
		return Graph_GML_grammar::parseGraph( GML_string );
    }
		
	
//##############################################################################

} // namespace ggl

