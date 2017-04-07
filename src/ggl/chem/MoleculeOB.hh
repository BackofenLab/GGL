#ifndef GGL_CHEM_MOLECULEOB_HH_
#define GGL_CHEM_MOLECULEOB_HH_


#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>


#include <sstream>
#include <iostream>
#include <limits>

#include "ggl/chem/Molecule.hh"


namespace ggl {
namespace chem {


	 /*! @brief OpenBabel molecule object port
	  *
	  * Wrapper class around an OpenBabel molecule object (OBMol).
	  *
	  * It utilizes the OpenBabel library (http://openbabel.org) for the
	  * conversion.
	  *
	  * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MoleculeOB
	{
	protected:

		//! current molecule working on
		OpenBabel::OBMol curMol;

	public:

		/*! Construction
		 *
		 * NOTE: the protons of the molecule are removed during construction.
		 *
		 * @param mol the molecule to be represented by this object
		 */
		MoleculeOB( const ggl::chem::Molecule & mol );

		/*! Access to the molar mass of the molecule
		 *
		 * @param implicitH flag whether hydrogens should be added implicitly by openbabel or not
		 *
		 * @return the standard molar mass given by IUPAC atomic masses (amu)
		 */
		double
		getMolWeight( const bool implicitH = true );

		/*! Access to the complete energy of the molecule, calculated via force field methods
		 *
		 * In error case (e.g. force field not loadable) NaN is returned.
		 *
		 * @return the complete energy of this molecule (in kJ/mol)
		 */
		double
		getMolEnergy();


		 /**
		  * Converts a Molecule into an OpenBabel molecule instance.
		  *
		  * NOTE: conversion is done via CML encoding and thus currently
		  * inefficient.
		  *
		  * @param mol the molecule to convert
		  * @return the according OBMol instance
		  */
		static
		OpenBabel::OBMol
		convert( const Molecule & mol );


		 /**
		  * Converts an OpenBabel molecule object into a Molecule instance.
		  *
		  * NOTE: only atom label and formal charge are considered for each
		  * atom.
		  *
		  * @param mol the molecule to convert
		  * @return the according Molecule instance
		  */
		static
		Molecule
		convert( const OpenBabel::OBMol & mol );

	};

} // namespace chem
} // namespace ggl

 // include function bodies
#include "ggl/chem/MoleculeOB.icc"


#endif /*GGL_CHEM_MOLECULEOB_HH_*/
