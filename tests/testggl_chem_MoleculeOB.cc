
#include <iostream>
#include <vector>
#include <string>

#include "ggl/chem/SMILESwriter.hh"
#include "ggl/chem/SMILESparser.hh"

#include <ggl/chem/MoleculeOB.hh>

using namespace ggl::chem;

void
molWeightTest (void) {

  std::vector< std::string > smilesInput;   // SMILES to test the getMolWeight() method with
  std::vector< double > lowerBoundWeightsWithHs;  // list of weights the molecules should have at least, to pass the test
  std::vector< double > lowerBoundWeightsWithoutHs;

  smilesInput.push_back("OCC1OC(O)C(O)C(O)C1O");
  lowerBoundWeightsWithHs.push_back(180.0);
  lowerBoundWeightsWithoutHs.push_back(165.0);

  smilesInput.push_back("[H]C(=O)[H]");
  lowerBoundWeightsWithHs.push_back(30.0);
  lowerBoundWeightsWithoutHs.push_back(27.0);

  smilesInput.push_back("[OH][CH3][CH2](O)[CH2](O)[CH2](O)[CH2](O)[CH2]O");
  lowerBoundWeightsWithHs.push_back(180.0);
  lowerBoundWeightsWithoutHs.push_back(165.0);

  for (size_t i = 0; i < smilesInput.size(); i++)
  {
    std::string smi = smilesInput[i];
    ggl::chem::Molecule molecule = SMILESparser::parseSMILES(smi).first;

    std::cout << "\ncheck '" << smi << "' from Molecule:" << "\n";

    MoleculeOB obmol2(molecule);
    std::cout << " weight():      "; std::cout.flush();
    double weight = obmol2.getMolWeight();
	std::cout << (weight > lowerBoundWeightsWithHs[i] ? "true" : "false") << "\n";
    std::cout << " weight(true):  "; std::cout.flush();
    weight = obmol2.getMolWeight(true);
	std::cout << (weight > lowerBoundWeightsWithHs[i] ? "true" : "false") << "\n";
    std::cout << " weight(false): "; std::cout.flush();
    weight = obmol2.getMolWeight(false);
	std::cout << (weight > lowerBoundWeightsWithoutHs[i] ? "true" : "false");
  }

  std::cout << std::endl;
}

void
molEnergyTest (void) {

  std::vector< std::string > smilesInput;   // SMILES to test the getMolWeight() method with
  std::vector< double > lowerBoundEnergyWithHs;  // list of weights the molecules should have at least, to pass the test

  smilesInput.push_back("OCC1OC(O)C(O)C(O)C1O");
  lowerBoundEnergyWithHs.push_back(1000.0);

  smilesInput.push_back("[H]C(=O)[H]");
  lowerBoundEnergyWithHs.push_back(0.0);

  smilesInput.push_back("[OH][CH3][CH2](O)[CH2](O)[CH2](O)[CH2](O)[CH2]O");
  lowerBoundEnergyWithHs.push_back(1000.0);

  for (size_t i = 0; i < smilesInput.size(); i++)
  {
    std::string smi = smilesInput[i];
    ggl::chem::Molecule molecule = SMILESparser::parseSMILES(smi).first;

    std::cout << "\ncheck '" << smi << "' from Molecule:" << "\n";

    MoleculeOB obmol2(molecule);
    std::cout << " energy(): "; std::cout.flush();
    double energy = obmol2.getMolEnergy();
	std::cout << (energy > lowerBoundEnergyWithHs[i] ? "true" : "false");
  }

  std::cout << std::endl;
}


void
convertTest (void) {

	  std::vector< std::string > smilesInput;   // SMILES to test the getMolWeight() method with
	  std::vector< double > lowerBoundEnergyWithHs;  // list of weights the molecules should have at least, to pass the test

	  smilesInput.push_back("OCC1OC(O)C(O)C(O)C1O");
	  smilesInput.push_back("[H]C(=O)[H]");
	  smilesInput.push_back("[OH][CH3][CH2](O)[CH2](O)[CH2](O)[CH2](O)[CH2]O");

	  for (size_t i=0; i<smilesInput.size(); ++i) {

		    ggl::chem::Molecule molecule = SMILESparser::parseSMILES(smilesInput.at(i)).first;
		    OpenBabel::OBMol molOB = MoleculeOB::convert( molecule );
		    ggl::chem::Molecule mol2 = MoleculeOB::convert( molOB );

		    std::cout <<"\n\tconvert from : '" <<SMILESwriter::getSMILES( molecule ) <<"'\n";
		    std::cout <<"\t          to : '" <<SMILESwriter::getSMILES( mol2 ) <<"'\n";
	  }

}

int main (int argc, char** argv) {

  molWeightTest();

  molEnergyTest();

  convertTest();

	return 0;
}
