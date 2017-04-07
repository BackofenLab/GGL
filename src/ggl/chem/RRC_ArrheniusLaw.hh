
#ifndef GGL_CHEM_RRC_ARRHENIUSLAW_HH_
#define GGL_CHEM_RRC_ARRHENIUSLAW_HH_

#include "ggl/chem/EnergyCalculation.hh"
#include "ggl/chem/ReactionRateCalculation.hh"

namespace ggl {
namespace chem {

///////////////////////////////////////////////////////////////////////////

	  /*! @brief Arrhenius law based reaction rate calculation
	   *
	   * Computes a reaction rate based on the Arrhenius law. To this end, it
	   * evaluates the energy difference of the educts compared to the products
	   * and uses the Boltzmann weight to predict a reaction rate.
	   *
	   * It uses the formula
	   *
	   *   rate k = A * exp( -deltaE / kT )
	   *
	   * using the Arrhenius constant (A), the energy difference of educts to
	   * products (deltaE), and the generalized temperature (kT).
	   *
	   * Note kT has to be set according to the used energy calculator. For
	   * its determination for a given temperature T, you can use the
	   * predefined constants within ggl::chem::EnergyCalculationConstants like
	   * the Boltzmann constant k_B or the gas constant R.
	   */
	class RRC_ArrheniusLaw : public ggl::chem::ReactionRateCalculation {

	protected:

		  //! the energy calculator to be used
		const EnergyCalculation * energyCalculation;

		  //! the generalized temperature to be applied within the Arrhenius law
		const double kT;

		  //! the Arrhenius constant to be applied within the Arrhenius law
		const double Arrhenius_constant_A;

	public:

		  /*!
		   * Construction
		   *
		   * @param energyCalculation the energy calculator to use to estimate
		   *        the energy of educts and products
		   * @param kT the generalized temperature to be applied within the
		   *        Arrhenius law
		   * @param Arrhenius_constant_A the Arrhenius constant to be applied
		   *        within the Arrhenius law
		   */
		RRC_ArrheniusLaw( const EnergyCalculation & energyCalculation,
						const double kT,
						const double Arrhenius_constant_A = 1.0 );


		  //! destruction
		virtual ~RRC_ArrheniusLaw();


		  /*!
		   * Calculates the reaction rate for a given Reaction based on the
		   * Arrhenius law.
		   *
		   * @param reaction the Reaction object to calculate the rate for
		   * @return the according reaction rate
		   */
		virtual
		double
		getRate( const Reaction & reaction ) const;

		  /*!
		   * The Arrhenius law calculation does not need the transition state
		   * information.
		   *
		   * @return false
		   */
		virtual
		bool
		needTransitionState( void ) const;

	};

///////////////////////////////////////////////////////////////////////////

} // namespace chem
}  // namespace ggl

#endif /* GGL_CHEM_RRC_ARRHENIUSLAW_HH_ */
