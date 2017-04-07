
#ifndef GGL_CHEM_EC_MOLECULEDECOMPOSITION_HH_
#define GGL_CHEM_EC_MOLECULEDECOMPOSITION_HH_

#include "ggl/chem/EnergyCalculation.hh"
#include "ggl/chem/MoleculeDecomposition.hh"

namespace ggl {
namespace chem {



	  /*! @brief Molecule energy ala Jankowski et al.
	   *
	   * An energy calculator implementation based on the molecule decomposition
	   * approach introduced by Jankowski et al. (2008)
	   * (see ggl::chem::MoleculeDecomposition)
	   */
	class EC_MoleculeDecomposition : public EnergyCalculation
	{

	protected:

		  //! the graph isomorphism matcher needed for the decomposition
		sgm::GraphMatching * fullMatcher;
		  //! the subgraph mono-morphism matcher needed for the decomposition
		sgm::SubGraphMatching * subMatcher;

		  //! the molecule decomposition handler used
		mutable MoleculeDecomposition * decompositionHandler;

	public:

		  //! construction
		EC_MoleculeDecomposition();

		  //! copy construction
		  //! @param toCopy the instance to make this a copy of
		EC_MoleculeDecomposition( const EC_MoleculeDecomposition& toCopy );

		  //! destruction
		virtual
		~EC_MoleculeDecomposition();


		  //! assignment operator
		  //! @param toCopy the instance to make this a copy of
		  //! @return the changed *this object
		virtual
		EC_MoleculeDecomposition &
		operator=( const EC_MoleculeDecomposition& toCopy );

		  /*!
		   * Computes the energy estimate in "kcal/mol" of a given molecule
		   * based of a decomposition of the molecule graph into components
		   * and a summation over all component contributions.
		   *
		   * @param mol the molecule to be analyzed
		   * @return an estimate of its energy in kcal/mol
		   */
		virtual
		double
		getEnergy( const Molecule & mol ) const;

	};


}  // namespace chem
} // namespace ggl

  // include member implementations
#include "ggl/chem/EC_MoleculeDecomposition.icc"

#endif /* GGL_CHEM_EC_MOLECULEDECOMPOSITION_HH_ */
