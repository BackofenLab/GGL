
#ifndef GGL_CHEM_ENERGYCALCULATIONCONSTANTS_HH_
#define GGL_CHEM_ENERGYCALCULATIONCONSTANTS_HH_

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/EnergyCalculationConstants.hh"

namespace ggl {
namespace chem {



	  /*! @brief Interface molecule energy calculation
	   *
	   * General interface for the prediction of a molecule's energy
	   * in kcal/mol units.
	   *
	   */
	class EnergyCalculation {

	public:

		  //! construction
		EnergyCalculation()
		{}

		  //! destruction
		virtual ~EnergyCalculation()
		{}


		  /*!
		   * Abstract function that computes the energy estimate in "kcal/mol"
		   * of a given molecule.
		   *
		   * @param mol the molecule to be analyzed
		   * @return an estimate of its energy in kcal/mol
		   */
		virtual
		double
		getEnergy( const Molecule & mol ) const = 0;


	};

}  // namespace chem
} // namespace ggl

#endif /* GGL_CHEM_ENERGYCALCULATIONCONSTANTS_HH_ */
