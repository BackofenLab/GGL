

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <exception>

#include "sgm/Graph_boost.hh"

#include "ggl/Graph.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/AP_NSPDK.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILESwriter.hh"

#include "dataTargetGraph_1.icc"

#include "utilPrintGraph_Interface.icc"



using namespace std;
using namespace ggl;
using namespace ggl::chem;

void testFeatures( AP_NSPDK & ap, const std::string molSMILES ) {

	cout <<" given     : '" <<molSMILES <<"'\n"; std::cout.flush();
	 // get molecule graph
	std::pair< Molecule, int > parse = SMILESparser::parseSMILES( molSMILES );
	cout <<" parse     : "<<parse.second <<std::endl; std::cout.flush();

	 // fill protons
	Molecule mol;
	mol = Molecule();
	MoleculeUtil::copy( parse.first, mol );
	MoleculeUtil::fillProtons( mol );
	cout <<" + protons : '"<<SMILESwriter::getSMILES( mol, false )<<"'" <<std::endl; std::cout.flush();

	cout <<" ->getFeatures( mol )" <<std::endl; std::cout.flush();
	std::vector< std::pair< sgm::RingReporter::RingList, nspdk::SVector> > features = ap.getFeatures( mol );

	  // create a dummy graph
	sgm::Graph_boost< ggl::chem::Molecule \
					, ggl::PropNodeLabel \
					, ggl::PropEdgeLabel \
					, ggl::PropNodeIndex \
					> molGraph(mol);

	for (size_t i=0; i<features.size(); i++) {
		std::cout <<" ring '";
		sgm::RingReporter::RingList::const_iterator cur = features.at(i).first.begin();
		sgm::RingReporter::RingList::const_iterator last = features.at(i).first.begin();
		  // initialize with first node label
		std::cout << molGraph.getNodeLabel( *cur );
		for (cur++; cur!=features.at(i).first.end();cur++,last++) {
			  // add bond
			std::cout << molGraph.getEdgesBegin( *last, *cur )->getEdgeLabel()
			  // ad next node label
					<< molGraph.getNodeLabel( *cur );
		}
		std::cout <<"' = " <<features.at(i).second <<std::endl;
	}
	std::cout <<std::endl; std::cout.flush();
}

void testMol( AP_NSPDK & ap, const std::string molSMILES ) {

	cout <<" given     : '" <<molSMILES <<"'\n"; std::cout.flush();
	 // get molecule graph
	Molecule mol, cmol;
	std::pair< Molecule, int > parse = SMILESparser::parseSMILES( molSMILES );
	cout <<" parse     : "<<parse.second <<std::endl; std::cout.flush();

	 // fill protons
	MoleculeUtil::copy( parse.first, mol );
	MoleculeUtil::fillProtons( mol );
	cmol = Molecule();
	MoleculeUtil::copy( mol, cmol );
	MoleculeUtil::compressHnodes( cmol );
	cout <<" + protons : '"<<SMILESwriter::getSMILES( cmol, false )<<"'" <<std::endl; std::cout.flush();

	try {
		 // rewrite molecule graph according to ring perception
		ap.correctAromaticity( mol, true );
	} catch (std::exception& e) {
		cout << " relabeling exception caught: " << e.what() << endl;
		cout << " according graph:\n" <<mol <<std::endl; std::cout.flush();
	}
//	  // create a dummy graph
//	sgm::Graph_boost< ggl::chem::Molecule \
//					, ggl::PropNodeLabel \
//					, ggl::PropEdgeLabel \
//					, ggl::PropNodeIndex \
//					> molGraph(parse.first);
//	cout <<" corrected graph :\n" << molGraph <<std::endl;

	 // get SMILES version of corrected molecule
	MoleculeUtil::compressHnodes( mol );
	cout <<" corrected : '"<<SMILESwriter::getSMILES( mol, false )<<"'" <<std::endl; std::cout.flush();
	cout <<" - protons : '"<<SMILESwriter::getSMILES( mol, true )<<"'\n"<<std::endl; std::cout.flush();

}


int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"         ggl::chem::AP_NSPDK  \n"
				<<"==============================================\n" 
				<<std::endl;

	std::vector< std::string > molNonArom;

	  // TODO more molecules
	molNonArom.push_back("O=C1C=CC(=O)C=C1");
	molNonArom.push_back("C1=CC=C1");
	molNonArom.push_back("NC(=O)CCCCC1CCSS1");
	molNonArom.push_back("O=C1CCCCC1");
	molNonArom.push_back("CC(C=CC=C(C)C=CC1=C(C)CCCC1(C)C)=CCO");
	molNonArom.push_back("CC(C=CC=C(C)C=CC1=C(C)CCCC1(C)C)=CC=O");
	molNonArom.push_back("[O-]C(=O)CCC1NC=NC1=O");
	molNonArom.push_back("CN1CC(=O)NC1=N");
	molNonArom.push_back("OC1CCCCC1");
//molNonArom.push_back("O=C1NC(=O)C2NC(=O)C(=O)NC=2N1"); // wrong in OpenBabel
//molNonArom.push_back("C1=CC=CC=CC=C1"); // wrong in OpenBabel
	molNonArom.push_back("C1CCCC#CCC1");

	std::vector< std::string > molArom;

	  // TODO more molecules
	molArom.push_back("c1c2c(nc[nH]2)ncn1");
//molArom.push_back("c1cc2C=Cc3cccc(c1)c23"); // wrong : OpenBabel + Marvin:general
	molArom.push_back("Oc1ccccc1");
	molArom.push_back("Nc1ncnc2[nH]cnc12");
	molArom.push_back("NC(=O)c1cccnc1");
	molArom.push_back("[O-]C(=O)C(=O)Cc1ccccc1");
	molArom.push_back("[O-]C(=O)C(=O)Cc1c[nH]c2ccccc12");
	molArom.push_back("[O-]C(=O)C(=O)Cc1cnc2ccccc12");
	molArom.push_back("COc1cc(cc(OC)c1O)C=CC([O-])=O");
	molArom.push_back("COc1cc(C=CCO)ccc1O");
	molArom.push_back("[H]C(=C([H])c1c[nH]cn1)C([O-])=O");
	molArom.push_back("CC(=O)n1ccnc1");
	molArom.push_back("COc1cc(C=CC=O)ccc1O");
	molArom.push_back("Cc1ncsc1CCO");
	molArom.push_back("Oc1ccc(cc1)C=CC(=O)c1ccc(O)cc1O");
	molArom.push_back("c1ccncc1");
	molArom.push_back("c1ccccc1-c2ccccc2");
	molArom.push_back("c1ccccc1c2ccccc2");


	molArom.push_back("c1[nH]c2c(n1)c(=O)[nH]c(n2)N");	// guanine
	molArom.push_back("c1nc2c(=O)[nH]c(nc2n1[CH]3[CH]([CH]([CH](O3)CO)O)O)N");	// guanosine
	molArom.push_back("c1nc2c(n1[CH]3C[CH]([CH](O3)CO)O)[nH]c(nc2=O)N");	// deoxyguanosine
	molArom.push_back("n1c(c2c(nc1)ncn2)N");	// adenine
	molArom.push_back("Nc1ncnc2n(cnc12)[CH]1O[CH](CO)[CH](O)[CH]1O");	// adenosine
	molArom.push_back("n2c1c(ncnc1n(c2)[CH]3O[CH]([CH](O)C3)CO)N");	// deoxyadenosine


	// tricky cases including sulphur (sulphur ring is predicted non-aromatic)
	molArom.push_back("Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)(O)=O)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(CC(O)=O)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(CCOP(O)(O)=O)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(C[CH]=O)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(CCOP(O)([O-])=O)c2C)c(N)n1");
	molArom.push_back("Cc1ncc(C[n+]2csc(CCOP(O)(=O)OP(O)([O-])=O)c2C)c(N)n1");

	molArom.push_back("Oc1ccc(cc1)c2cn3c(nc(Cc4ccccc4)c3=O)c(Cc5ccccc5)[nH]2");
	molArom.push_back("CC[CH](C)c1nc2c(CCCNC(N)=N)[nH]c(cn2c1=O)c3c[nH]c4ccccc43");
	molArom.push_back("OS(=O)(=O)Oc1ccc(Cc2nc3c(Cc4ccccc4)[nH]c(cn3c2=O)c5ccc(OS(O)(=O)=O)cc5)cc1");
	molArom.push_back("Oc1ccc(Cc2nc3c(Cc4ccccc4)[nH]c(cn3c2=O)c5ccc(O)cc5)cc1");

	// weird cases with biiiig aromatic ring containing several smaller ones
	// --> currently not working since covering large aromatic ring is not predicted and ignored, thus backup handling does not work
//	molArom.push_back("Cc1c(CCC(O)=O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)c(C)c4C=C)c(C)c3CCC(O)=O");
//	molArom.push_back("Cc1c(CCC(O)=O)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)c(C)c5CCC(O)=O)c(CCC(O)=O)c4C)c(CCC(O)=O)c3C");
//	molArom.push_back("Cc1c(CCC(O)=O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CCC(O)=O)c5C)c(CCC(O)=O)c4C)c(CCC(O)=O)c3C");


	const std::set < std::string > & models = AP_NSPDK_Model::getAvailableModels();

	for ( std::set < std::string >::const_iterator modelID = models.begin(); modelID != models.end(); ++modelID) {

		std::cout <<"\n################################################\n\n";
		std::cout <<"load model '"<<*modelID<<"'\n"; std::cout.flush();
		AP_NSPDK apSVM(*modelID);

		testFeatures( apSVM, "c1cc2C=Cc3cccc(c1)c23" );

		std::cout <<"\n**********************************\n\n";
		std::cout <<"\n testing non-aromatic molecules :\n"<<std::endl;
		std::vector<std::string>::const_iterator mol;
		for ( mol=molNonArom.begin(); mol != molNonArom.end(); mol++ ) {
			try {
				testMol( apSVM, *mol );
			} catch (std::exception & ex) {
				std::cout <<"  --> aromaticity perception failed : "
						<<ex.what() <<std::endl;
			}
		}

		std::cout <<"\n**********************************\n\n";
		std::cout <<"\n testing aromatic molecules :\n"<<std::endl;
		for ( mol=molArom.begin(); mol != molArom.end(); mol++ ) {
			try {
				testMol( apSVM, *mol );
			} catch (std::exception & ex) {
				std::cout <<"  --> aromaticity perception failed : "
						<<ex.what() <<std::endl;
			}
		}
	}



	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
