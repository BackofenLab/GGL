
#ifndef GGL_CHEM_SMILES_PARSER_HH_
#define GGL_CHEM_SMILES_PARSER_HH_

#include <utility>
#include <string>

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeUtil.hh"
#include "ggl/chem/ChemRule.hh"

namespace ggl {
  namespace chem {


	  /*! @brief SMILES molecule parser
	   *
	   *  This class is a wrapper for the ggl::chem::SMILESparser that
	   *  parses SMILES via their BNF grammar. See there for further details.
	   *
	   *  @author Martin Mann (c) 2012 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class SMILESparser
	{
	public:
	    
	    
	      /**
	       *  Parses a SMILES string and generates a graph representation of the
	       *  molecule
	       *  @param SMILES_string the string to parse
	       *  @return pair.first = the graph encoding of the molecule
	       *          pair.second = -1 if parsing was successful,
	       *          in error case it returns the string position that caused
	       *          the parsing error
		   *  @throw std::invalid_argument in case a check fails
		   */
	    static
	    std::pair< Molecule, int >
	    parseSMILES( const std::string & SMILES_string )
	    	throw (std::invalid_argument);
	    

	      /**
	       *  Parses a SMILES string and generates a graph representation of the
	       *  molecule
	       *  @param SMILES_string the string to parse
	       *  @param groups a container that holds group IDs where
	       *        each matching node has to be replaced by the according
	       *        mapped subgraph
	       *  @return pair.first = the graph encoding of the molecule
	       *          pair.second = -1 if parsing was successfull,
	       *          in error case it returns the string position that caused
	       *          the parsing error
		   *  @throw std::invalid_argument in case a check fails
		   */
	    static
	    std::pair< Molecule, int >
	    parseSMILES( const std::string & SMILES_string, const GroupMap & groups )
	    	throw (std::invalid_argument);


	      /**
	       *  Parses a SMILES string of multiple molecules separated by '.'
	       *  and generates a graph representation of the molecules
	       *  @param SMILES_string the string to parse
	       *  @return vector of according molecule graphs
		   *  @throw std::invalid_argument in case a parsing fails
		   */
	    static
	    std::vector< Molecule * >
	    parseMultiSMILES( const std::string & SMILES_string )
	    	throw (std::invalid_argument);


	      /**
	       *  Parses a SMILES string of multiple molecules separated by '.'
	       *  and generates a graph representation of the molecules
	       *  @param SMILES_string the string to parse
	       *  @param groups a container that holds group IDs where
	       *        each matching node has to be replaced by the according
	       *        mapped subgraph
	       *  @return vector of according molecule graphs
		   *  @throw std::invalid_argument in case a parsing fails
		   */
	    static
	    std::vector< Molecule * >
	    parseMultiSMILES( const std::string & SMILES_string, const GroupMap & groups )
	    	throw (std::invalid_argument);


	      /**
	       * Parses a reaction SMILES string and uses the atom mapping encoded
	       *  in the atom label class information to generate a graph grammar
	       *  rule core graph, thus, each atom in the educts has to
	       *  have a unique class ID and a corresponding atom (with same ID) in
	       *  the products.
	       *
	       *  An example '[C:13][C-2:2].[O:3]>>[O:3]=[C:2][C:13]'
	       *
	       *  For protons part of complex labels like [CH3:1] are automatically
	       *  given class IDs while it is expected that the adjacent atom (with
	       *  given class ID) is present both in educt and product and features
	       *  in both the same number of adjacent protons given without class
	       *  ID within the complex atom label. Otherwise, an exception is
	       *  raised.
	       *
	       *  NOTE: the reaction is not checked for sanity or chemical
	       *  correctness. Furthermore, no constraints etc. are generated.
	       *
	       *  @param SMILES_string the string to parse
	       *  @param pruneClassID whether or not the classID information is to
	       *         be pruned from the atom labels
	       *  @return the parsed rule core graph
		   *  @throw std::invalid_argument with according description in case
		   *         a check fails
		   */
	    static
	    ChemRule::CoreGraph
	    parseReactionSMILES( const std::string & SMILES_string
	    					, const bool pruneClassID = true )
	    	throw (std::invalid_argument);

	protected :


	      /**
	       * Prunes all aromatic bonds that are not part of a ring.
	       *
	       * @param mol the molecule to correct
	       */
	    static
	    void
	    pruneAromaticNonRingBonds( Molecule& mol );

  }; // class SMILESparser
  

  
  } // namespace chem
} // namespace ggl



#endif /*GGL_CHEM_SMILES_PARSER_HH_*/
