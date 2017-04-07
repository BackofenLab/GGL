#ifndef GGL_RULE_GML_PARSER_HH_
#define GGL_RULE_GML_PARSER_HH_

#include <utility>

#include "ggl/Rule.hh"
#include "ggl/Rule_GML_error.hh"

namespace ggl {



	/*! @brief Graph grammar rule parser
	 *
	 * Wrapper for ggl::Rule_GML_grammar BNF grammar parser that parses a
	 * GML string representation of a ggl::Rule object. See there for further
	 * details.
	 * 
	 * @author Martin Mann (c) 2012 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */ 
	class Rule_GMLparser
	{
		
	public:
	    
		  /*!
	       * Parses a GML string and generates a Rule object
	       *
	       * @param GML_string the string to parse
	       *
	       * @return pair.first = the graph encoding of the molecule
	       *         pair.second = -1 if parsing was successful,
	       *         in error case it returns the string position that caused
	       *         the parsing error
	       *
	       * @throw ggl::Rule_GML_error if parsing error occur
		   */
	    static
	    std::pair< Rule, int >
	    parseRule( const std::string & GML_string ) throw (Rule_GML_error);
	    
		  /*!
	       * Parses a GRAPH GML string that encodes a compacted RULE GML
	       * and generates a Rule object.
	       *
	       * NOTE: rule constraints, its label, or wildcard information are not
	       * parsed, just the nodes and edges!
	       *
	       * @param GML_string the string to parse
	       *
	       * @return pair.first = the graph encoding of the molecule
	       *         pair.second = -1 if parsing was successful,
	       *         in error case it returns the string position that caused
	       *         the parsing error
	       *
	       * @throw ggl::Rule_GML_error if parsing error occur
		   */
	    static
	    std::pair< Rule, int >
	    parseCompactRule( const std::string & GML_string ) throw (Rule_GML_error);

	};


} // namespace ggl

#endif /*RULE_GML_PARSER_HH_*/
