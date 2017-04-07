
#include "ggl/chem/RRC_ArrheniusLaw.hh"
#include "ggl/chem/SMILESparser.hh"

namespace ggl {
namespace chem {


///////////////////////////////////////////////////////////////////////////////


	RRC_ArrheniusLaw
	::RRC_ArrheniusLaw( const EnergyCalculation & energyCalculation,
						const double kT,
						const double Arrhenius_constant_A )
	  :	energyCalculation(&energyCalculation)
		, kT(kT)
		, Arrhenius_constant_A(Arrhenius_constant_A)
	{
	}


///////////////////////////////////////////////////////////////////////////////


	RRC_ArrheniusLaw
	::~RRC_ArrheniusLaw()
	{
	}


///////////////////////////////////////////////////////////////////////////////


	double
	RRC_ArrheniusLaw
	::getRate( const Reaction & reaction ) const
	{
		double rate = 0.0;

		  // get energy of educts
		double energyEducts = 0.0;
		for (Reaction::Metabolite_Container::const_iterator m=reaction.metabolites.begin();
				m != reaction.metabolites.end(); ++m )
		{
			  // get molecule graph
			std::pair< Molecule, int > mol = SMILESparser::parseSMILES( *m );
			assert(mol.second == -1);
			  // fill protons
			MoleculeUtil::fillProtons( mol.first );
			  // compute energy and sum up
			energyEducts += energyCalculation->getEnergy( mol.first );
		}

		  // get energy of products
		double energyProducts = 0.0;
		for (Reaction::Product_Container::const_iterator p=reaction.products.begin();
				p != reaction.products.end(); ++p )
		{
			  // get molecule graph
			std::pair< Molecule, int > mol = SMILESparser::parseSMILES( *p );
			assert(mol.second == -1);
			  // fill protons
			MoleculeUtil::fillProtons( mol.first );
			  // compute energy and sum up
			energyProducts += energyCalculation->getEnergy( mol.first );

		}

		  // calculate rate based on Arrhenius law
		rate = Arrhenius_constant_A * exp( -(energyProducts-energyEducts) / kT );

		return rate;
	}


///////////////////////////////////////////////////////////////////////////////


	bool
	RRC_ArrheniusLaw
	::needTransitionState( void ) const
	{
		  // no transition state needed
		return false;
	}


///////////////////////////////////////////////////////////////////////////////


}  // namespace chem
} // namespace ggl
