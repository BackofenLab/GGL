
#include "ggl/chem/MoleculeComponent_GMLparser.hh"
#include "ggl/chem/MoleculeComponent_GML_grammar.hh"

namespace ggl {
namespace chem {

//##############################################################################

	std::pair< MoleculeComponent, int >
	MoleculeComponent_GMLparser
    ::parseGML( const std::string & GML_string ) throw (std::invalid_argument)
    {
		// forward call
		return MoleculeComponent_GML_grammar::parseGML( GML_string );
    }
		
	
//##############################################################################

} // namespace chem
} // namespace ggl

