#ifndef GGL_CHEM_GS_MOLCHECK_HH_
#define GGL_CHEM_GS_MOLCHECK_HH_

#include <iostream>

#include "ggl/chem/GS_chem.hh"
#include "ggl/chem/Molecule.hh"
#include "ggl/chem/AromaticityPerception.hh"

namespace ggl {
namespace chem {


	////////////////////////////////////////////////////////////////////////////

	  /*! @brief Molecule correction graph storage
	   *
	   * A graph storage that checks each reported molecule graph to consistency
	   * and does needed corrections.
	   *
	   * Applied corrections are e.g.
	   *  - aromaticity perception and relabeling
	   *
	   * The corrected molecule is than reported to the given successive graph
	   * storage.
	   *
	   * @author Martin Mann (c) 2011 http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */
	class GS_MolCheck : public GS_chem {

	protected:

		 //! the graph storage to which checked molecules are forwarded to
		GS_chem* nextGS;

		 //! aromaticity perception handler to relabel aromatic and non-aromatic
		 //! rings
		AromaticityPerception * aromaticityPerception;

		 //! whether or not exceptions are to be ignored or forwarded
		bool ignoreExceptions;

	public:
		
		  /*!
		   * Construction.
		   *
		   * @param gs the graph storage to report the checked and relabeled
		   *          molecules to
		   * @param aromaticityPerception the aromaticity perception instance
		   *          to be used to perform the aromaticity relabeling;
		   * @param ignoreExceptions if true: exceptions raised during molecule
		   *          check are ignored; otherwise the exception is forwarded
		   */
		GS_MolCheck( GS_chem & gs
					, const AromaticityPerception & aromaticityPerception
					, const bool ignoreExceptions = true );

		  /*!
		   * Copy construction
		   * @param toCopy the instance to make this object a copy of
		   */
		GS_MolCheck( const GS_MolCheck & toCopy );


		  /*!
		   * Assignment operator
		   * @param toCopy the instance to make this object a copy of
		   * @return the altered object *this
		   */
		GS_MolCheck &
		operator=( const GS_MolCheck & toCopy );


		  //! destruction
		virtual ~GS_MolCheck();
		
		  /*!
		   * Checks a given molecule graph for inconsistent labeling and
		   * corrects it accordingly. The corrected molecule is than forwarded
		   * to the successive graph storage.
		   *
		   * Note: errors during the aromaticity rewrites etc. will also result
		   * in the disposal of the molecule without further notification!
		   *
		   * @param mol the molecule graph object to check and add.
		   *
		   * @throw std::runtime_error in case the molecule correction fails
		   *          and ignoreExceptions is set to false; otherwise no
		   *          exception is thrown
		   */
		virtual
		void
		addMolecule( const Molecule & mol );
	};
	
	////////////////////////////////////////////////////////////////////////////


}} // namespace ggl

#include "ggl/chem/GS_MolCheck.icc"

#endif /*GGL_CHEM_GS_MOLCHECK_HH_*/
