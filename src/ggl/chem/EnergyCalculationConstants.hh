
#ifndef GGL_CHEM_ENERGYCALCULATION_HH_
#define GGL_CHEM_ENERGYCALCULATION_HH_

#include "ggl/chem/Molecule.hh"


namespace ggl {
namespace chem {


	  /*! @brief Constants for energy calculations
	   *
	   */
	class EnergyCalculationConstants {

	public:

		  //! The gas constant R given in unit "kcal/(mol*K)"
		static const double Gas_constant_R;

		  //! The Avogardo constant N_A given in unit "1/mol"
		static const double Avogardo_constant_NA;

		  //! The Boltzmann constant k_B = R / N_A, i.e. gas constant (R)
		  //! divided by the Avogardo constant (N_A), given in unit "kcal/K"
		static const double Boltzmann_constant_kB;

	};

}  // namespace chem
} // namespace ggl

#endif /* GGL_CHEM_ENERGYCALCULATION_HH_ */
