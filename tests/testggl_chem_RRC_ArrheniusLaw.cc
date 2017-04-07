
#include <iostream>
#include <sstream>

#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILESwriter.hh"
#include "ggl/chem/MoleculeUtil.hh"
#include "ggl/chem/EC_MoleculeDecomposition.hh"
#include "ggl/chem/RRC_ArrheniusLaw.hh"

#include "sgm/Graph_boost.hh"

#include "utilPrintGraph_Interface.icc"

using namespace ggl;
using namespace ggl::chem;


const std::string SMILES[] = {

//		"ENERGY","SMILES", // "ENTRY","NAMES","GROUP_DECOMPOSITION"

		"-480.93","OP([O-])(=O)OP([O-])([O-])=O", // 8,"diphosphate|pyrophosphate|ppi|","pyrophosphate : 1"
		"-18.97","[H]N([H])[H]", // 9,"ammonia|nh4|","NH4 : 1"
		"-112.69","CC(=O)C([O-])=O", // 15,"pyruvate|pyr|",">C=O : 1 | -COO : 1 | -CH3 : 1 | Origin : 1 | OCCO : 1"
		"-188.9","[O-]C(=O)CCC(=O)C([O-])=O", // 17,"alpha-ketoglutarate|akg|",">C=O : 1 | -CH2- : 2 | -COO : 2 | Origin : 1 | OCCO : 1"
		"-32.05","[H]OO[H]", // 18,"hydrogen peroxide|","H2O2 : 1"
		"-190.52","OC(=O)CC(=O)C([O-])=O", // 20,"oxaloacetate|oaa|",">C=O : 1 | -CH2- : 1 | -COO : 2 | Origin : 1 | OCCO : 1"
		"-87.73","[N+]CC([O-])=O", // 21,"glycine|gly|","-COO : 1 | -CH2- : 1 | -NH3 : 1 | Origin : 1"
		"-87.73","[NH3+]CC([O-1])=O", // 21,"glycine|gly|","-COO : 1 | -CH2- : 1 | -NH3 : 1 | Origin : 1"
		"-162.96","[O-]C(=O)CCC([O-])=O", // 23,"succinate|succ|","-COO : 2 | -CH2- : 2 | Origin : 1"

	"" // marks the end of the array
};

std::string getSMILES(const std::string & smiles, EnergyCalculation & energyCalc, double& energyAdd );

int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"      ggl::chem::RRC_ArrheniusLaw             \n"
				<<"==============================================\n" 
				<<std::endl;

	/*!
	 * - create reaction
	 * - print energy sums
	 * - set temperature
	 * - calculate rate
	 */

	const double kT = ggl::chem::EnergyCalculationConstants::Gas_constant_R * 298.15;
	  // setup energy calculator
	EC_MoleculeDecomposition energyCalc;
	  // setup reaction rate calculator
	RRC_ArrheniusLaw rateCalculator( energyCalc, kT );

	std::cout <<"\n used temperature in K = " <<(kT/ggl::chem::EnergyCalculationConstants::Gas_constant_R) <<"\n" <<std::endl;

	  // setup dummy reaction
	Reaction reaction;
	reaction.rule_id = "dummy reaction";
	double metaboliteEnergies = 0.0, productEnergies = 0.0;
	std::string smiles;
	Molecule mol, mol2;

	//////////////////////////////////////////////////////////////////////////

	reaction.rule_id = "ATP + H2O -> ADP + phosphate";

	  // setup metabolites
	reaction.metabolites.clear(); metaboliteEnergies = 0.0;

	reaction.metabolites.insert( getSMILES("Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C1O", energyCalc, metaboliteEnergies ) );
	reaction.metabolites.insert( getSMILES("[H]O[H]", energyCalc, metaboliteEnergies ) );

	  // setup products
	reaction.products.clear(); productEnergies = 0.0;
	reaction.products.insert( getSMILES("Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C1O", energyCalc, productEnergies ) );
	reaction.products.insert( getSMILES("OP(O)(O)=O", energyCalc, productEnergies ) );

	std::cout <<"\n next reaction  : " <<reaction <<std::endl;


	std::cout <<" E(metabolites) = "<< metaboliteEnergies <<" kcal/mol"<<std::endl;
	std::cout <<" E(products)    = "<< productEnergies <<" kcal/mol" <<std::endl;
	std::cout <<" delta E        = "<< (productEnergies-metaboliteEnergies) <<" kcal/mol" <<std::endl;
	std::cout <<" rate(reaction) = "<< rateCalculator.getRate( reaction ) <<std::endl;

	//////////////////////////////////////////////////////////////////////////

//	  // setup metabolites
//	reaction.metabolites.clear();metaboliteEnergies = 0.0;
//	smiles="[H]N([H])[H]"; reaction.metabolites.insert(smiles);metaboliteEnergies += getEnergy(energyCalc, smiles);
//	smiles="[H]OO[H]"; reaction.metabolites.insert(smiles);metaboliteEnergies += getEnergy(energyCalc, smiles);
//
//	  // setup products
//	reaction.products.clear();productEnergies = 0.0;
//	smiles="[NH3+]CC([O-1])=O"; reaction.products.insert(smiles);productEnergies += getEnergy(energyCalc, smiles);
//
//	std::cout <<"\n next reaction  : " <<reaction <<std::endl;
//
//
//	std::cout <<" E(metabolites) = "<< metaboliteEnergies <<std::endl;
//	std::cout <<" E(products)    = "<< productEnergies <<std::endl;
//	std::cout <<" rate(reaction) = "<< rateCalculator.getRate( reaction ) <<std::endl;

	//////////////////////////////////////////////////////////////////////////

	std::cout	<<"\n"
				<<"===============  END TEST  ===================\n"
				<<std::endl;
	return 0;
}


std::string getSMILES(const std::string & smiles, EnergyCalculation & energyCalc, double& energyAdd )
{
	Molecule mol=SMILESparser::parseSMILES(smiles).first;
	MoleculeUtil::fillProtons(mol);
	energyAdd += energyCalc.getEnergy(mol);
	MoleculeUtil::compressHnodes(mol);
	std::string newSmiles = SMILESwriter::getSMILES(mol);

	return newSmiles;

}

