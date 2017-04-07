
#include "ggl/chem/GS_MolCheck.hh"


namespace ggl {
namespace chem {

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


	GS_MolCheck&
	GS_MolCheck::
	operator=( const GS_MolCheck & toCopy )
	{
		  // copy graph storage
		this->nextGS = toCopy.nextGS;
		  // copy aromaticity perception handler
		if (this->aromaticityPerception != NULL)
			delete this->aromaticityPerception;
		this->aromaticityPerception = toCopy.aromaticityPerception->clone();
		  // copy exception handling
		this->ignoreExceptions = toCopy.ignoreExceptions;
		  // return this
		return *this;
	}

////////////////////////////////////////////////////////////////////////////



	void
	GS_MolCheck::
	addMolecule( const Molecule & graph )
	{

		try {
			Molecule mol(graph);
			  // try aromaticity rewrite
			aromaticityPerception->correctAromaticity( mol, false );

			  // forward call to successive graph storage
			nextGS->add(mol);

		} catch (std::exception& ex) {
			if (!ignoreExceptions) {
				throw std::runtime_error(ex.what());
			}
			// do nothing : most important : dont store buggy molecule
			return;
		}
	}

////////////////////////////////////////////////////////////////////////////

}} //namespace
