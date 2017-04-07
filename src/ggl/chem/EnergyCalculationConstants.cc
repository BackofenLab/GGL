
#include "ggl/chem/EnergyCalculationConstants.hh"

#include <cmath>

namespace ggl {
namespace chem {


	  // The gas constant R given in unit "kcal/(mol*K)"
	const double EnergyCalculationConstants::Gas_constant_R = 8.314462175 * 2.3901 * pow(10,-4);

	  // The Avogardo constant N_A given in unit "1/mol"
	const double EnergyCalculationConstants::Avogardo_constant_NA = 6.0221412927 * pow(10,23);

	  // The Boltzmann constant k_B = R / N_A, i.e. gas constant (R)
	  // divided by the Avogardo constant (N_A), given in unit "kcal/K"
	const double EnergyCalculationConstants::Boltzmann_constant_kB = (8.314462175 * 2.3901 / 6.0221412927) * pow(10,-19);


}  // namespace chem
} // namespace ggl

