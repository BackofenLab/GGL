#ifndef GGL_CHEM_SMILESWRITEROB_HH_
#define GGL_CHEM_SMILESWRITEROB_HH_


#include <sstream>
#include <iostream>
#include <memory>

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeUtil.hh"

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>



namespace ggl {
 namespace chem {

	 /*! @brief Molecule to SMILES writer via OpenBabel
	  *
	  * Utility class to generate a canonical SMILES string from a molecule
	  * graph representation. It expects atom (node) and edge (bond) label
	  * following the Daylight's SMILES description. See 
	  * ggl::chem::SMILES_grammar for further details.
	  * 
	  * It utilizes the OpenBabel library (http://openbabel.org) for the 
	  * conversion.
	  * 
	  * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	  * 
	  */
	class SMILESwriterOB
	{
	
	protected:

		  //! the OpenBabel conversion object used
		static std::auto_ptr<OpenBabel::OBConversion> converter;
	
	public:
		SMILESwriterOB();
		
		  /*! Generates a canonical SMILES string of the given graph 
		   * representation of a molecule.
		   * 
		   * NOTE: THE FUNCTIONALITY IS INCOMPLETE, i.e. only a few atom and
		   * bond types are possible! See class description!
		   * 
		   * @param mol the molecule graph to parse
		   * @return a canonical SMILES string representing the given molecule
		   */
		static
		std::string
		getSMILES( const Molecule& mol)
		{
			  // setup converter if not done yet
			if (converter.get() == 0) {
				  // create converter
				converter.reset( new OpenBabel::OBConversion(NULL,NULL));
				  // setup conversion scheme
				converter->SetInAndOutFormats("CML","CAN");
			}

			  // get CML representation of the molecule
			std::stringstream molCML, outSMILES;
			MoleculeUtil::convertCML( mol, molCML );

			  // convert from CML to canonical SMILES
			converter->Convert( &molCML, &outSMILES );

			  // trim whitespaces
			std::string out = outSMILES.str();
			size_t from = out.find_first_not_of(" \t\n");
				assert(from != out.npos);
			size_t to = out.find_last_not_of(" \t\n");
				assert(to != out.npos);
			  // return conversion result
			return out.substr(from, to-from+1);
		}
		
		
	};
	
	  // create member
	std::auto_ptr<OpenBabel::OBConversion>
	SMILESwriterOB::converter = std::auto_ptr<OpenBabel::OBConversion>();

 } // namespace chem
} // namespace ggl


#endif /*GGL_CHEM_SMILESWRITEROB_HH_*/
