
#include <iostream>
#include <vector>

#include <sgm/Graph_boost.hh>
#include <ggl/chem/MoleculeUtil.hh>
#include <ggl/chem/SMILESparser.hh>

void
labelTest (void) {
	
	using namespace ggl::chem;
	
	std::vector< std::string > validAtomLabel;
	
	  // add all supported atom label
	const MoleculeUtil::AtomDataMap & lm = MoleculeUtil::getAtomData();
	for ( MoleculeUtil::AtomDataMap::const_iterator i=lm.begin(); i!=lm.end(); i++) {
		validAtomLabel.push_back(i->first);
	}
	  // add special label
	
	  // H
	validAtomLabel.push_back("CH");
	validAtomLabel.push_back("CH2");
	  // charge
	validAtomLabel.push_back("H+");
	validAtomLabel.push_back("O-2");
	validAtomLabel.push_back("O-10");
	validAtomLabel.push_back("OH-");
	validAtomLabel.push_back("OH-2");
	validAtomLabel.push_back("CH3-");
	validAtomLabel.push_back("CH2-2");
	  // class
	validAtomLabel.push_back("H:1");
	validAtomLabel.push_back("O-:2");
	validAtomLabel.push_back("O-2:2");
	validAtomLabel.push_back("OH:21");
	validAtomLabel.push_back("OH2:213");
	validAtomLabel.push_back("OH-:21");
	validAtomLabel.push_back("CH3-:2");
	validAtomLabel.push_back("CH2-2:12");
	
	std::cout <<"\n check SYNTACTICALLY VALID atom label :\n\n";
	for (size_t i=0; i< validAtomLabel.size(); i++) {
		std::cout	<<" check '"<<validAtomLabel.at(i) <<"' : "; std::cout.flush();
		std::cout	<< (MoleculeUtil::isValidAtomLabel( validAtomLabel.at(i) )?"true":"false"); std::cout.flush();
		std::cout	<<"\n  atom     = " <<MoleculeUtil::getAtom( validAtomLabel.at(i) ); std::cout.flush();
		std::cout	<<"\n  protons  = " <<MoleculeUtil::getProtons( validAtomLabel.at(i) ); std::cout.flush();
		std::cout	<<"\n  charge   = " <<MoleculeUtil::getCharge( validAtomLabel.at(i) ); std::cout.flush();
		std::cout	<<"\n  class    = " <<MoleculeUtil::getClass( validAtomLabel.at(i) ); std::cout.flush();
		std::cout	<<std::endl;
		const MoleculeUtil::AtomLabelData* ad = MoleculeUtil::getAtomData( validAtomLabel.at(i) );
		std::cout	<<"  atom data = ";
		if (ad == NULL)
			std::cout <<"NULL";
		else {
			std::cout	<<" valence = " <<(int)ad->valence 
						<<" atomic number = " <<(int)ad->atomicNumber
						<<" aromatic ? " <<(int(ad->isAromatic)==1?"true":"false");
		}
		std::cout	<<std::endl;

		std::cout	<<validAtomLabel.at(i)
					<<" == "
					<<MoleculeUtil::getComplexAtomLabel(
								MoleculeUtil::getAtom( validAtomLabel.at(i) )
								, MoleculeUtil::getProtons( validAtomLabel.at(i) )
								, MoleculeUtil::getCharge( validAtomLabel.at(i) )
								, MoleculeUtil::getClass( validAtomLabel.at(i) )
								, false)
					<<" == "
					<<MoleculeUtil::getComplexAtomLabel(
								MoleculeUtil::getAtom( validAtomLabel.at(i) )
								, MoleculeUtil::getProtons( validAtomLabel.at(i) )
								, MoleculeUtil::getCharge( validAtomLabel.at(i) )
								, MoleculeUtil::getClass( validAtomLabel.at(i) )
								, true)
					<<"\n";
		// TODO test MoleculeUtil::getComplexAtomLabel(..) via string comparison
	}
	
	  // creating invalid labels
	std::vector< std::string > invalidAtomLabel;
	
	  // white spaces
	invalidAtomLabel.push_back("");
	invalidAtomLabel.push_back(" ");
	invalidAtomLabel.push_back(" H");
	invalidAtomLabel.push_back("H ");
	invalidAtomLabel.push_back("C a");
	invalidAtomLabel.push_back("[C ]");
	invalidAtomLabel.push_back("[C] ");
	invalidAtomLabel.push_back(" [C]");
	invalidAtomLabel.push_back(" [ C]");
	invalidAtomLabel.push_back(" [ C ]");
	  // SMILES brackets left
	invalidAtomLabel.push_back("[C]");
	invalidAtomLabel.push_back("[CH]");
	invalidAtomLabel.push_back("[C+]");
	invalidAtomLabel.push_back("[CH+]");
	invalidAtomLabel.push_back("[C:3]");
	invalidAtomLabel.push_back("[CH:3]");
	invalidAtomLabel.push_back("[C+:3]");
	invalidAtomLabel.push_back("[CH+:3]");
	invalidAtomLabel.push_back("[C");
	invalidAtomLabel.push_back("C]");
	  // protons
	invalidAtomLabel.push_back("HH");
	invalidAtomLabel.push_back("H2");
	invalidAtomLabel.push_back("CHH");
	invalidAtomLabel.push_back("CH12");
	  // charge
	invalidAtomLabel.push_back("C+H");
	invalidAtomLabel.push_back("C-2H");
	invalidAtomLabel.push_back("C-2H2");
	invalidAtomLabel.push_back("C-2+");
	invalidAtomLabel.push_back("C+-");
	invalidAtomLabel.push_back("C+123");
	  // class
	invalidAtomLabel.push_back("H:");
	invalidAtomLabel.push_back("O:-");
	invalidAtomLabel.push_back("O:2-");
	invalidAtomLabel.push_back("O:2-2");
	invalidAtomLabel.push_back("O:21H");
	invalidAtomLabel.push_back("O:21H2");
	invalidAtomLabel.push_back("OH:21-");
	invalidAtomLabel.push_back("O:21H-");
	invalidAtomLabel.push_back("C:2H3-");
	invalidAtomLabel.push_back("CH3:2-");
	invalidAtomLabel.push_back("C:12H2-2");
	invalidAtomLabel.push_back("CH2:12-2");
	  // strange things
	invalidAtomLabel.push_back("CC");
	invalidAtomLabel.push_back("CHC");
	invalidAtomLabel.push_back("CH2C");
	invalidAtomLabel.push_back("C+C");
	invalidAtomLabel.push_back("C+2C");
	invalidAtomLabel.push_back("CH+2C");
	invalidAtomLabel.push_back("CH2+2C");
	invalidAtomLabel.push_back("C:2C");
	invalidAtomLabel.push_back("CH:2C");
	invalidAtomLabel.push_back("CH2:2C");
	invalidAtomLabel.push_back("CH2+:2C");
	invalidAtomLabel.push_back("CH2+2:2C");
	
	
	std::cout <<"\n check SYNTACTICALLY INVALID atom label :\n\n";
	for (size_t i=0; i< invalidAtomLabel.size(); i++) {
		std::cout	<<" check '"<<invalidAtomLabel.at(i) <<"' : "
					<< (MoleculeUtil::isValidAtomLabel( invalidAtomLabel.at(i) )?"true":"false")
					<<std::endl;
	}
}


void
fillProtonTest()
{
	using namespace ggl;
	using namespace ggl::chem;

	std::vector< std::string > molecules;
	molecules.push_back("[H+]");
	molecules.push_back("HH");
	molecules.push_back("[HH]");
	molecules.push_back("O");
	molecules.push_back("[OH-]");
	molecules.push_back("[O-]");
	molecules.push_back("[O+]");
	molecules.push_back("[OH3+]");
	molecules.push_back("O=C1C=CC(=O)C=C1");
	molecules.push_back("CC[H:1]");
	  // coordinative bonds
	molecules.push_back("C^C");
//	molecules.push_back("[Fe+2](^O)^O");

	size_t consistencyCode = 0;
	for (size_t i=0; i<molecules.size(); i++) {
		std::cout <<"\n fillProtons('" <<molecules.at(i) <<"') :" <<std::endl;

		std::pair< Molecule, int > parse = SMILESparser::parseSMILES( molecules.at(i) );
		std::cout <<" parsing result : " <<parse.second <<std::endl;

		std::cout <<" before:\n" <<sgm::Graph_boost<Molecule>(parse.first) <<std::endl;
		std::cout <<" isConsistent( .. ) = "; std::cout.flush();
		consistencyCode = MoleculeUtil::isConsistent(parse.first);
		std::cout <<(consistencyCode == MoleculeUtil::C_Consistent?"true":"false") <<std::endl;
		if (consistencyCode != MoleculeUtil::C_Consistent) {
			MoleculeUtil::decodeConsistencyStatus( consistencyCode, std::cout );
		}

		MoleculeUtil::fillProtons( parse.first );
		std::cout <<" after:\n" <<sgm::Graph_boost<Molecule>(parse.first) <<std::endl;
		std::cout <<" isConsistent( .. ) = "; std::cout.flush();
		consistencyCode = MoleculeUtil::isConsistent(parse.first);
		std::cout <<(consistencyCode == MoleculeUtil::C_Consistent?"true":"false") <<std::endl;
		if (consistencyCode != MoleculeUtil::C_Consistent) {
			MoleculeUtil::decodeConsistencyStatus( consistencyCode, std::cout );
		}

		std::cout <<"\n compressHnodes( ... ) :" <<std::endl;
		MoleculeUtil::compressHnodes( parse.first );
		std::cout <<" after:\n" <<sgm::Graph_boost<Molecule>(parse.first) <<std::endl;

	}

}


int main (int argc, char** argv) {
	
	labelTest();
	fillProtonTest();
	
	return 0;
}
