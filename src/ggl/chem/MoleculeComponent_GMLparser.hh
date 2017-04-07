#ifndef GGL_CHEM_MOLECULECOMPONENT_GML_PARSER_HH_
#define GGL_CHEM_MOLECULECOMPONENT_GML_PARSER_HH_

#include <utility>
#include <string>
#include <stdexcept>

#include "ggl/chem/MoleculeComponent.hh"

namespace ggl {
namespace chem {


	/*! @brief MoleculeComponent parser
	 *
	 * This class is a wrapper for the ggl::chem::MoleculeComponent_GML_grammar
	 * BNF grammar parser. It parses a GML string representation of a
	 * ggl::MoleculeDecomposition::MoleculeComponent object. This includes its
	 * properties as well as the additional constraints needed for matching.
	 * See ggl::chem::MoleculeComponent_GML_grammar for further details
	 * 
	 * @author Martin Mann (c) 2012 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */ 
	class MoleculeComponent_GMLparser
	{
		
	public:
	    
	      //! Parses a GML string and generates a MoleculeComponent::PatternGraph object
	      //! @param GML_string the string to parse
	      //! @return pair.first = the graph encoding of the molecule
	      //!         pair.second = -1 if parsing was successfull,
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
	      //! @throw std::invalid_argument in case a check fails
	    static
	    std::pair< MoleculeComponent, int >
	    parseGML( const std::string & GML_string ) throw (std::invalid_argument);
	};


} // namespace chem
} // namespace ggl


#endif /*GGL_CHEM_MOLECULECOMPONENT_GML_PARSER_HH_*/
