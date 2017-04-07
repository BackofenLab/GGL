
#include <iostream>
#include <sstream>

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeDecomposition.hh"
#include "ggl/chem/SMILESparser.hh"

#include "sgm/Graph_boost.hh"
#include "sgm/GM_vf2.hh"
#include "sgm/SGM_vf2.hh"

#include "utilPrintGraph_Interface.icc"

using namespace ggl;
using namespace ggl::chem;

// SMILES REFERENCE : http://metacyc.org/

const std::string SMILES[] = {

//		"ENERGY","SMILES", // "ENTRY","NAMES","GROUP_DECOMPOSITION"
		"-480.93","OP([O-])(=O)OP([O-])([O-])=O", // 8,"diphosphate|pyrophosphate|ppi|","pyrophosphate : 1"
		"-18.97","[H]N([H])[H]", // 9,"ammonia|nh4|","NH4 : 1"
		"-18.97","[NH4+]", // 9,"ammonia(aq)|nh4|","NH4 : 1"
		"-112.69","CC(=O)C([O-])=O", // 15,"pyruvate|pyr|",">C=O : 1 | -COO : 1 | -CH3 : 1 | Origin : 1 | OCCO : 1"
		"-188.9","[O-]C(=O)CCC(=O)C([O-])=O", // 17,"alpha-ketoglutarate|akg|",">C=O : 1 | -CH2- : 2 | -COO : 2 | Origin : 1 | OCCO : 1"
		"-32.05","[H]OO[H]", // 18,"hydrogen peroxide|","H2O2 : 1"
		"-190.52","OC(=O)CC(=O)C([O-])=O", // 20,"oxaloacetate|oaa|",">C=O : 1 | -CH2- : 1 | -COO : 2 | Origin : 1 | OCCO : 1"
		"-87.73","[N+]CC([O-])=O", // 21,"glycine|gly|","-COO : 1 | -CH2- : 1 | -NH3 : 1 | Origin : 1"
		"-87.73","[NH3+]CC([O-1])=O", // 21,"glycine|gly|","-COO : 1 | -CH2- : 1 | -NH3 : 1 | Origin : 1"
		"-162.96","[O-]C(=O)CCC([O-])=O", // 23,"succinate|succ|","-COO : 2 | -CH2- : 2 | Origin : 1"
		"-111.04","[H]C(=O)C([O-])=O", // 26,"glyoxylate|glx|","-COO : 1 | -CH=O : 1 | Origin : 1 | OCCO : 1"
		"-83.9","[H]C([O-])=O", // 32,"formate|for|","formate : 1"
		"-177.97","[O-]S([O-])(=O)=O", // 33,"sulfate|so4|","sulfate : 1"
		"-31.2","[H]C([H])=O", // 38,"formaldehyde|fald|","formaldehyde : 1"
		"-316.08","[O-]P([O-])(=O)OC(=C)C([O-])=O", // 40,"phosphoenolpyruvate|pep|","-OPO3 : 1 | =C< : 1 | -COO : 1 | =CH2 : 1 | Origin : 1 | OCCC : 1"
		"-33.4","[H]C(C)=O", // 48,"acetaldehyde|acald|","acetaldehyde : 1"
		"-48.7","NC(N)=O", // 50,"urea|","urea : 1"
		"-86.11","[N+]CCC([O-])=O", // 54,"beta-alanine|","-CH2- : 2 | -COO : 1 | -NH3 : 1 | Origin : 1"
		"-61.07","O=C1NC=CC(=O)N1", // 58,"uracil|","ring >C=O : 2 | ring =CH- : 2 | amide : 2 | ring -NH- : 2 | Origin : 1 | OCCC : 1"
		"-111.07","CCC(=O)C([O-])=O", // 59,"2-oxobutanoate|",">C=O : 1 | -COO : 1 | -CH2- : 1 | -CH3 : 1 | Origin : 1 | OCCO : 1"
		"12.19","C[N+](C)(C)CCO", // 61,"choline|","-CH2- : 2 | >N< : 1 | -CH3 : 3 | -OH : 1 | Origin : 1"
		"-116.18","OCC(O)CO", // 62,"glycerol|glyc|",">CH- : 1 | -CH2- : 2 | -OH : 3 | Origin : 1"
		"-84.87","CC(C)CC([N+])C([O-])=O", // 66,"leucine|leu|",">CH- : 2 | -CH2- : 1 | -COO : 1 | -NH3 : 1 | -CH3 : 2 | Origin : 1"
		"-439.84","CC(=C)CCOP([O-1])(=O)OP([O-])([O-1])=O", // 68,"isopentenyl-diphosphate|","-OPO2- : 1 | -OPO3 : 1 | -CH2- : 2 | =C< : 1 | -CH3 : 1 | =CH2 : 1 | Origin : 1"
		"-634.01","Nc1ncnc2n(cnc12)C1CC(O)C(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)O1", // 70,"dATP","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 3 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring -CH2- : 1 | ring =C< : 1 | -NH2 : 1 | -CH2- : 1 | -OH : 1 | -OPO2- : 2 | -OPO3 : 1 | Origin : 1 | heteroaromatic : 2"
		"-41.9","CO", // 71,"methanol|","methanol : 1"
		"-87.92","CC([N+])C([O-])=O", // 72,"alanine|ala|",">CH- : 1 | -COO : 1 | -CH3 : 1 | -NH3 : 1 | Origin : 1"
		"-9.9","Oc1ccccc1", // 76,"phenol|","aromatic ring =C< : 1 | aromatic ring =CH- : 5 | -OH : 1 | Origin : 1"
		"78.85","Nc1ncnc2[nH]cnc12", // 77,"adenine|","two fused rings >C= : 2 | ring =C< : 1 | ring =N- : 3 | ring -NH- : 1 | -NH2 : 1 | ring =CH- : 2 | Origin : 1 | heteroaromatic : 2"
		"5.55","NC(=O)c1cccnc1", // 80,"nicotinamide|","ring =C< : 1 | ring =CH- : 4 | >C=O : 1 | ring =N- : 1 | amide : 1 | -NH2 : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 1"
		"-280.44","OC(CC([O-])=O)(CC([O-])=O)C([O-])=O", // 83,"citrate|cit|",">C< : 1 | -CH2- : 2 | -COO : 3 | -OH : 1 | Origin : 1"
		"-122.98","OCC([O-])=O", // 84,"glycolate|glyclt|","-COO : 1 | -CH2- : 1 | -OH : 1 | Origin : 1"
		"-85.13","CCC([O-])=O", // 85,"propanoate|ppa|","-COO : 1 | -CH2- : 1 | -CH3 : 1 | Origin : 1"
		"-75.82","[O-]C(=O)C(=O)Cc1ccccc1", // 86,"keto-phenylpyruvate|phenylpyruvate|","aromatic ring =C< : 1 | -CH2- : 1 | aromatic ring =CH- : 5 | >C=O : 1 | -COO : 1 | Origin : 1 | OCCO : 1"
		"-148.92","OCC(=O)C([O-])=O", // 87,"3-hydroxypyruvate|",">C=O : 1 | -COO : 1 | -CH2- : 1 | -OH : 1 | Origin : 1 | OCCO : 1"
		"-108.16","OCC(=O)CO", // 90,"dihydroxyacetone|dha|glycerone",">C=O : 1 | -CH2- : 2 | -OH : 2 | Origin : 1"
		"-38.9","NO", // 92,"hydroxylamine|","hydroxylamine : 1"
		"-38.52","CC(C)=O", // 95,"acetone|","acetone : 1"
		"-111.88","[H]C(=O)CC([O-])=O", // 98,"3-oxopropanoate|","-CH2- : 1 | -COO : 1 | -CH=O : 1 | Origin : 1"
		"-301.65","CC(=O)OP([O-])([O-])=O", // 99,"acetyl-phosphate|acetyl-p|","-CO-OPO3- : 1 | -CH3 : 1 | Origin : 1"
		"-110.26","[O-]C(=O)CCC=O", // 100,"4-oxobutanoate|","-CH2- : 2 | -COO : 1 | -CH=O : 1 | Origin : 1"
		"-439.18","CC(C)=CCOP([O-])(=O)OP([O-])([O-])=O", // 102,"dimethylallyl-diphosphate|","-OPO2- : 1 | -OPO3 : 1 | -CH2- : 1 | =CH- : 1 | =C< : 1 | -CH3 : 2 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"28.4","Nc1nc2[nH]cnc2c(=O)[nH]1", // 105,"guanine|","two fused rings >C= : 2 | ring =N- : 2 | ring -NH- : 2 | ring >C=O : 1 | ring =C< : 1 | ring =CH- : 1 | -NH2 : 1 | amide : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-21.54","NC(=O)CCCCC1CCSS1", // 106,"lipoamide|","ring -CH< : 1 | ring -CH2- : 2 | ring -S- : 2 | -CH2- : 4 | >C=O : 1 | amide : 1 | -NH2 : 1 | Origin : 1"
		"-64.07","CCCCCCCCCCCCCCCC([O-])=O", // 107,"palmitate|","-CH2- : 14 | -COO : 1 | -CH3 : 1 | Origin : 1"
		"-123.17","CC(O)C([O-])=O", // 109,"lactate|lac-L|",">CH- : 1 | -COO : 1 | -CH3 : 1 | -OH : 1 | Origin : 1"
		"-0.35","[H]C(=O)c1ccccc1", // 112,"benzaldehyde|","aromatic ring =C< : 1 | aromatic ring =CH- : 5 | -CH=O : 1 | Origin : 1 | OCCC : 1"

// only 1 heteroaromatic ring by jankowski
//	"23.12","O=c1[nH]cnc2nc[nH]c12", // 113,"hypoxanthine|Hypoxanthine|hypoxanthine(0,4)|hypoxanthine|hxan|","two fused rings >C= : 2 | ring >C=O : 1 | ring =N- : 2 | ring -NH- : 2 | ring =CH- : 2 | amide : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-70.28","[H]C(=O)CO", // 115,"glycolaldehyde|","-CH2- : 1 | -CH=O : 1 | -OH : 1 | Origin : 1"
		"-142.48","C([O-])(=O)C1(=CC(=O)NC(N1)=O)", // 122,"orotate|","ring =C< : 1 | ring =CH- : 1 | -COO : 1 | ring >C=O : 2 | amide : 2 | ring -NH- : 2 | Origin : 1 | OCCC : 2"
		"-65.59","CN(CC([O-])=O)C(N)=[N+]", // 125,"creatine|Creatine|creatine(0,9)|creatine|","-N< : 1 | -CH2- : 1 | =C< : 1 | -CH3 : 1 | -COO : 1 | -NH2 : 1 | =NH2 : 1 | Origin : 1"
		"-279.02","OC(C(CC([O-])=O)C([O-])=O)C([O-])=O", // 127,"isocitrate|icit|",">CH- : 2 | -CH2- : 1 | -COO : 3 | -OH : 1 | Origin : 1"
		"-49.77","[O-]C(=O)C(=O)Cc1c[nH]c2ccccc12", // 130,"indole-3-pyruvate|","aromatic ring fused to nonaromatic ring >C=  : 2 | ring =C< : 1 | aromatic ring =CH- : 4 | -CH2- : 1 | ring =CH- : 1 | ring -NH- : 1 | >C=O : 1 | -COO : 1 | Origin : 1 | OCCO : 1 | heteroaromatic : 1"
		"-112.07","CC(C=O)C([O-])=O", // 137,"2-methyl-3-oxopropanoate|",">CH- : 1 | -COO : 1 | -CH=O : 1 | -CH3 : 1 | Origin : 1"
		"-257.01","[N+]CCOP([O-])([O-])=O", // 145,"O-phosphorylethanolamine|","-OPO3 : 1 | -CH2- : 2 | -NH3 : 1 | Origin : 1"
		"-90.93","C12(NC(=O)NC=1C(=O)NC(=O)N2)", // 146,"urate|CHEBI:17775","two fused rings >C= : 2 | ring >C=O : 3 | amide : 4 | ring -NH- : 4 | Origin : 1 | OCCC : 1"
		"-14.2","O=C1CCCCC1", // 150,"cyclohexanone|","ring >C=O : 1 | ring -CH2- : 5 | Origin : 1"
		"-32.27","[O-]C(=O)C=Cc1ccccc1", // 152,"trans-cinnamate|","aromatic ring =C< : 1 | =CH- : 2 | aromatic ring =CH- : 5 | -COO : 1 | Origin : 1 | OCCC : 1 | CCCC : 1"
		"70.9","CC(C=CC=C(C)C=CC1=C(C)CCCC1(C)C)=CCO", // 160,"retinol|","ring =C< : 2 | ring >C< : 1 | =CH- : 6 | ring -CH2- : 3 | -CH3 : 5 | =C< : 2 | -CH2- : 1 | -OH : 1 | Origin : 1 | CCCC : 4"
		"-121.41","COc1cc(cc(OC)c1O)C=CC([O-])=O", // 163,"sinapate|","aromatic ring =C< : 4 | aromatic ring =CH- : 2 | -O- : 2 | -OH : 1 | -CH3 : 2 | =CH- : 2 | -COO : 1 | Origin : 1 | OCCC : 1 | CCCC : 1"
		"-65.83","C[N+](C)(C)CC(O)CC([O-])=O", // 164,"carnitine|","-CH2- : 2 | >N< : 1 | >CH- : 1 | -CH3 : 3 | -OH : 1 | -COO : 1 | Origin : 1"
		"-166.94","NC(=O)NC(NC(N)=O)C([O-])=O", // 168,"allantoate|",">CH- : 1 | -COO : 1 | >C=O : 2 | amide : 4 | -NH- : 2 | -NH2 : 2 | Origin : 1"
		"-85.66","[N+]CCCC([N+])C([O-])=O", // 171,"ornithine|orn|",">CH- : 1 | -COO : 1 | -CH2- : 3 | -NH3 : 2 | Origin : 1"
		"-59.99","[H]C(=O)C(C)=O", // 173,"methylglyoxal|",">C=O : 1 | -CH=O : 1 | -CH3 : 1 | Origin : 1 | OCCO : 1"
		"-8.28","OCc1ccccc1", // 176,"benzyl-alcohol|","aromatic ring =C< : 1 | aromatic ring =CH- : 5 | -CH2- : 1 | -OH : 1 | Origin : 1"
		"-106.7","[H]C(=O)C(O)CO", // 180,"glyceraldehyde|glyald|",">CH- : 1 | -CH2- : 1 | -CH=O : 1 | -OH : 2 | Origin : 1"
		"-27.35","NC(=O)CCCCC(S)CCS", // 181,"dihydrolipoamide|","-CH2- : 6 | >CH- : 1 | -SH : 2 | >C=O : 1 | amide : 1 | -NH2 : 1 | Origin : 1"
		"1.47","CSC", // 182,"dimethyl sulfide|","-S- : 1 | -CH3 : 2 | Origin : 1"
		"-75.84","NC(=[N+])NCC([O-])=O", // 183,"guanidinoacetic acid|","-COO : 1 | -CH2- : 1 | -NH- : 1 | =C< : 1 | -NH2 : 1 | =NH2 : 1 | Origin : 1"
		"-200.31","C[N+](C)(C)CCOP([O-])([O-])=O", // 186,"choline phosphate|","-CH2- : 2 | >N< : 1 | -CH3 : 3 | -OPO3 : 1 | Origin : 1"
		"-51.81","COc1cc(C=CCO)ccc1O", // 187,"coniferol|","aromatic ring =C< : 3 | -OH : 2 | aromatic ring =CH- : 3 | -O- : 1 | =CH- : 2 | -CH3 : 1 | -CH2- : 1 | Origin : 1 | CCCC : 1"
		"-114.09","CC(=O)NCCCC([O-])=O", // 193,"N-acetyl-L-alpha-amino-n-butyrate|CHEBI:17645","-CH2- : 3 | -COO : 1 | >C=O : 1 | -CH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-63.98","[O-]C(=O)CCCCC1SSCC1", // 200,"alpha-lipoate|","ring -CH< : 1 | ring -CH2- : 2 | ring -S- : 2 | -CH2- : 4 | -COO : 1 | Origin : 1"
		"-124.15","[N+]C(CO)C([O-])=O", // 201,"serine|",">CH- : 1 | -COO : 1 | -CH2- : 1 | -NH3 : 1 | -OH : 1 | Origin : 1"
		"-26.85","[H]C(=C([H])c1c[nH]cn1)C([O-])=O", // 202,"urocanate|","ring =C< : 1 | =CH- : 2 | ring =CH- : 2 | ring =N- : 1 | ring -NH- : 1 | -COO : 1 | Origin : 1 | OCCC : 1 | CCCC : 1 | heteroaromatic : 1"
		"-72.12","CC(O)C(C)=O", // 204,"acetoin|",">CH- : 1 | >C=O : 1 | -CH3 : 2 | -OH : 1 | Origin : 1"

// wrong SMILES
//		"-34.72","CC(O)C(=O)C1(CNC2(=NC(N)=NC(=O)[C-]2(N=1)))", // 206,"sepiapterin|","two fused rings >C= : 2 | ring -NH- : 2 | ring =N- : 2 | ring >C=O : 1 | ring -CH2- : 1 | ring =C< : 2 | -NH2 : 1 | >C=O : 1 | >CH- : 1 | -CH3 : 1 | -OH : 1 | amide : 1 | Origin : 1 | OCCC : 1 | OCCN : 1 | CCNC : 2"

		"-20.76","OC1CCCCC1", // 208,"cyclohexanol|","ring -CH< : 1 | ring -CH2- : 5 | -OH : 1 | Origin : 1"
		"-145.2","CC(C([O-])=O)=C(C)C([O-])=O", // 216,"dimethylmaleate|","-COO : 2 | =C< : 2 | -CH3 : 2 | Origin : 1 | OCCC : 2"
		"-70.47","[H]C(=O)C(C)O", // 217,"lactaldehyde|",">CH- : 1 | -CH=O : 1 | -CH3 : 1 | -OH : 1 | Origin : 1"
		"-119.74","OCCCC([O-])=O", // 223,"4-hydroxybutyrate|","-CH2- : 3 | -COO : 1 | -OH : 1 | Origin : 1"
		"-121.36","[O-]C(=O)CCO", // 224,"3-hydroxypropanoate|b-Hydroxypropionate|hydroxypropionate(-1,5)|hydroxypropionate|b-hydroxypropionate|","-CH2- : 2 | -COO : 1 | -OH : 1 | Origin : 1"
		"-115.71","CC(=O)NCCC([O-])=O", // 226,"N-acetyl-L-alanine|CHEBI:16682","-COO : 1 | -CH2- : 2 | >C=O : 1 | -CH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-199.38","OC(CCC([O-])=O)C([O-])=O", // 228,"2-hydroxyglutarate|",">CH- : 1 | -CH2- : 2 | -COO : 2 | -OH : 1 | Origin : 1"
		"-226.94","OC(CC(=O)C([O-])=O)C([O-])=O", // 237,"4-hydroxy-2-oxoglutarate|",">CH- : 1 | -COO : 2 | -CH2- : 1 | -OH : 1 | >C=O : 1 | Origin : 1 | OCCO : 1"
		"-149.92","[O-]C(C=O)C([O-])=O", // 240,"2-hydroxy-3-oxopropanoate|",">CH- : 1 | -COO : 1 | -CH=O : 1 | -OH : 1 | Origin : 1"
		"-115.3","[O-]C(=O)C(=O)Cc1ccc(O)cc1", // 241,"4-hydroxyphenylpyruvate|",">C=O : 1 | -COO : 1 | -CH2- : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -OH : 1 | Origin : 1 | OCCO : 1"
		"-121.55","CC(CO)C([O-])=O", // 243,"3-hydroxy-2-methylpropanoate|",">CH- : 1 | -COO : 1 | -CH2- : 1 | -CH3 : 1 | -OH : 1 | Origin : 1"
		"-79.76","OCCO", // 251,"ethylene-glycol|","-CH2- : 2 | -OH : 2 | Origin : 1"
		"-143.7","[O-]C(=O)C=CC([O-])=O", // 252,"maleate|fumarate|fum|","-COO : 2 | =CH- : 2 | Origin : 1 | OCCC : 2"
		"-30.81","[H]C(=O)CCC", // 253,"butanal|butyraldehyde|","-CH2- : 2 | -CH=O : 1 | -CH3 : 1 | Origin : 1"
		"-16.11","CCCCCC=CCC=CCCCCCCCC([O-])=O", // 256,"9-cis,12-cis-octadecadienoate|","-CH2- : 12 | =CH- : 4 | -COO : 1 | -CH3 : 1 | Origin : 1"
		"22.2","N#C[S-]", // 265,"thiocyanate|","#C- : 1 | #N : 1 | -Sneg : 1 | Origin : 1"
		"-43.72","CC(C)O", // 269,"isopropanol|",">CH- : 1 | -CH3 : 2 | -OH : 1 | Origin : 1"
		"-114.87","NC(=[N+])NCCS([O-])(=O)=O", // 272,"taurocyamine|","-SO3 : 1 | -CH2- : 2 | -NH- : 1 | =C< : 1 | -NH2 : 1 | =NH2 : 1 | Origin : 1"
		"-144.58","OCC(O)C(=O)CO", // 274,"erythrulose|",">CH- : 1 | >C=O : 1 | -CH2- : 2 | -OH : 3 | Origin : 1"
		"78.83","CC(C=CC=C(C)C=CC1=C(C)CCCC1(C)C)=CC=O", // 275,"11-cis-retinal|all-trans-retinal|","ring =C< : 2 | ring >C< : 1 | =CH- : 6 | ring -CH2- : 3 | -CH3 : 5 | =C< : 2 | -CH=O : 1 | Origin : 1 | OCCC : 1 | CCCC : 4"
		"-125.76","CC([N+])(CO)C([O-])=O", // 276,"2-methylserine|",">C< : 1 | -COO : 1 | -CH2- : 1 | -CH3 : 1 | -NH3 : 1 | -OH : 1 | Origin : 1"
		"-113.53","CC(=O)CC([O-])=O", // 277,"3-oxobutanoate|acetoacetate|acac|b-ketobutyrate|","-CH2- : 1 | -COO : 1 | -CH3 : 1 | >C=O : 1 | Origin : 1"
		"-69.79","[O-]C(=O)CCCCC(S)CCS", // 278,"dihydro-alpha-lipoate|","-CH2- : 6 | >CH- : 1 | -SH : 2 | -COO : 1 | Origin : 1"
		"-143.75","CC(C([O-])=O)C(=C)C([O-])=O", // 282,"methylitaconate|",">CH- : 1 | =C< : 1 | -COO : 2 | -CH3 : 1 | =CH2 : 1 | Origin : 1 | OCCC : 1"
		"-302.57","[H]C(=O)OP(O)(O)=O", // 284,"formyl-phosphate|","formylphosphate : 1"
		"3.35","CC(=O)n1ccnc1", // 289,"N-acetylimidazole|",">C=O : 1 | ring =CH- : 3 | -CH3 : 1 | ring =N- : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | heteroaromatic : 1"
		"-57.81","C[N+](C)(C)CC(=O)CC([O-])=O", // 296,"3-dehydrocarnitine|",">N< : 1 | -CH2- : 2 | -CH3 : 3 | >C=O : 1 | -COO : 1 | Origin : 1"
		"-43.88","COc1cc(C=CC=O)ccc1O", // 299,"coniferyl-aldehyde|","aromatic ring =C< : 3 | -OH : 1 | aromatic ring =CH- : 3 | -O- : 1 | =CH- : 2 | -CH3 : 1 | -CH=O : 1 | Origin : 1 | OCCC : 1 | CCCC : 1"
		"-80.78","[O-]C(=O)CNC=[N+]", // 303,"N-formiminoglycine|","-CH2- : 1 | -COO : 1 | -NH- : 1 | =CH- : 1 | =NH2 : 1 | Origin : 1"
		"-70.87","OC(=Cc1ccccc1)C([O-])=O", // 304,"enol-phenylpyruvate|","aromatic ring =C< : 1 | =CH- : 1 | aromatic ring =CH- : 5 | =C< : 1 | -COO : 1 | -OH : 1 | Origin : 1 | OCCC : 1 | CCCC : 1"
		"-101.66","CCC(=O)OC(CC([O-])=O)C[N+](C)(C)C", // 310,"propionylcarnitine|",">CH- : 1 | -O-CO- : 1 | -CH2- : 3 | >N< : 1 | -COO : 1 | -CH3 : 4 | Origin : 1"
		"-39.83","C(=O)c1cccc(O)c1", // 314,"3-hydroxybenzaldehyde|","aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -CH=O : 1 | -OH : 1 | Origin : 1 | OCCC : 1"
		"-272.38","[O-]C(=O)CNC(=[N+])NP([O-])([O-])=O", // 318,"phosphoguanidinoacetate|","-OPO2 : 1 | -NH- : 2 | =C< : 1 | =NH2 : 1 | -CH2- : 1 | -COO : 1 | Origin : 1"
		"-118.57","O=C1NC(=O)C2NC(=O)C(=O)NC=2N1", // 319,"tetrahydroxypteridine|","two fused rings >C= : 2 | ring >C=O : 4 | amide : 4 | ring -NH- : 4 | Origin : 1 | OCCO : 1 | OCCC : 1"
		"-157.39","CC([N+]C(C)C([O-])=O)C([O-])=O", // 320,"2,2'-iminodipropanoate|",">CH- : 2 | -COO : 2 | -NH2- : 1 | -CH3 : 2 | Origin : 1"

// only one heteroaromatic ring by jankowski ... most probably wrong SMILES
//	"-110.61","[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1", // 322,"xanthine-8-carboxylate|","two fused rings >C= : 2 | ring -NH- : 3 | ring >C=O : 2 | ring =N- : 1 | ring =C< : 1 | -COO : 1 | amide : 2 | Origin : 1 | OCCC : 1 | OCCN : 1 | heteroaromatic : 1"

		"-47.76","OCc1cccc(O)c1", // 323,"3-hydroxybenzyl-alcohol|","aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -CH2- : 1 | -OH : 2 | Origin : 1"
		"-86.31","[O-]C(=O)CCC1NC=NC1=O", // 339,"4,5-dihydro-4-oxo-5-imidazolepropanoate|","ring -CH< : 1 | ring >C=O : 1 | ring -NH- : 1 | -CH2- : 2 | ring =N- : 1 | ring =CH- : 1 | -COO : 1 | Origin : 1 | OCNC : 1"
		"-127.36","[O-]C(=O)CNC(=O)OCc1ccccc1", // 341,"benzyloxycarbonylglycine|","aromatic ring =C< : 1 | -CH2- : 2 | aromatic ring =CH- : 5 | -O-CO- : 1 | -COO : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-65.26","C[S+](C)CCC([O-])=O", // 348,"dimethylpropiothetin|","-COO : 1 | -CH2- : 2 | -S< : 1 | -CH3 : 2 | Origin : 1"
		"-20.93","CCCCCCC=CC=CCCCCCCCC([O-])=O", // 350,"9-cis,11-trans-octadecadienoate|","-CH2- : 12 | =CH- : 4 | -COO : 1 | -CH3 : 1 | Origin : 1 | CCCC : 1"
		"-132.27","[O-]C(=O)CC1OC(=O)C=C1", // 352,"2,5-dihydro-5-oxofuran-2-acetate|","ring -CH< : 1 | ring -O-CO- : 1 | ring =CH- : 2 | -CH2- : 1 | -COO : 1 | Origin : 1 | OCCC : 1"
		"-7.11","Cc1ncsc1CCO", // 359,"4-methyl-5-2'-hydroxyethyl-thiazole|","ring =C< : 2 | -CH2- : 2 | ring -S- : 1 | ring =N- : 1 | -CH3 : 1 | ring =CH- : 1 | -OH : 1 | Origin : 1 | heteroaromatic : 1"
		"-187.28","[O-]C(=O)CCCC(=O)C([O-])=O", // 381,"2-oxoadipate|","-CH2- : 3 | >C=O : 1 | -COO : 2 | Origin : 1 | OCCO : 1"
		"-65.96","Oc1ccc(cc1)C=CC(=O)c1ccc(O)cc1O", // 395,"2',4,4'-trihydroxychalcone|","aromatic ring =C< : 5 | >C=O : 1 | aromatic ring =CH- : 7 | =CH- : 2 | -OH : 3 | Origin : 1 | OCCC : 2 | CCCC : 1"
		"-140.26","OC([O-])=O", // 399,"hco3|bicarbonate|","HCO3 : 1"
		"-32.43","CCC=O", // 400,"propanal|","-CH3 : 1 | -CH2- : 1 | -CH=O : 1 | Origin : 1"
		"-146.69","[O-]C(=O)C1=NC(=O)NC(=O)N1", // 406,"oxonic-acid|","ring >C=O : 2 | ring =C< : 1 | -COO : 1 | ring =N- : 1 | amide : 2 | ring -NH- : 2 | Origin : 1 | OCCN : 1 | OCNC : 1"
		"33.64","Nc1ccccc1", // 412,"aniline|","aromatic ring =CH- : 5 | aromatic ring =C< : 1 | -NH2 : 1 | Origin : 1"
		"-25.26","CC(=O)OCC[N+](C)(C)C", // 413,"acetylcholine|","-CH3 : 4 | >N< : 1 | -CH2- : 2 | -O-CO- : 1 | Origin : 1"
		"-41.91","CCCO", // 449,"n-propanol|","-CH3 : 1 | -CH2- : 2 | -OH : 1 | Origin : 1"
		"-311.41","OP([O-])(=O)NC(=[N+])NCCS([O-])(=O)=O", // 452,"phosphotaurocyamine|","-OPO2 : 1 | -NH- : 2 | =C< : 1 | =NH2 : 1 | -CH2- : 2 | -SO3 : 1 | Origin : 1"
		"-45.29","[N+]CC(N)=O", // 459,"glycinamide|","-NH3 : 1 | -CH2- : 1 | >C=O : 1 | amide : 1 | -NH2 : 1 | Origin : 1"
		"-84.04","[N+]CCCCC([N+])C([O-])=O", // 465,"lysine|","-COO : 1 | -CH2- : 4 | >CH- : 1 | -NH3 : 2 | Origin : 1"
		"-80.46","[O-]C(=O)CNC(=O)Cc1ccccc1", // 471,"phenylacetylglycine|","-COO : 1 | -CH2- : 2 | >C=O : 1 | aromatic ring =C< : 1 | aromatic ring =CH- : 5 | amide : 1 | -NH- : 1 | Origin : 1"
		"-320.66","OCC(=O)COP([O-])([O-])=O", // 476,"glycerone-phosphate|","-OH : 1 | -CH2- : 2 | >C=O : 1 | -OPO3 : 1 | Origin : 1"
		"-78.65","CCCCCCC([O-])=O", // 490,"heptanoate|","-CH3 : 1 | -CH2- : 5 | -COO : 1 | Origin : 1"
		"-23.3","[O-]C#N", // 523,"cyanate|","#C- : 1 | #N : 1 | -Oneg : 1 | Origin : 1"
		"-161","OC(=O)C([O-])=O", // 524,"oxalate|","oxalate : 1"
		"-4.07","CC", // 525,"ethane|","ethane : 1"
		"-88.29","CC([O-])=O", // 526,"acetate|ac|","acetate : 1"
		"-43.44","CCO", // 527,"ethanol|etoh|","ethanol : 1"
		"-83.51","CCCC([O-])=O", // 530,"butyrate|but|","-CH2- : 2 | -COO : 1 | -CH3 : 1 | Origin : 1"
		"-81.89","CCCCC([O-])=O", // 531,"valerate|","-CH2- : 3 | -COO : 1 | -CH3 : 1 | Origin : 1"
		"-64.98","[O-]C(=O)C=C", // 532,"acrylate|","-COO : 1 | =CH- : 1 | =CH2 : 1 | Origin : 1 | OCCC : 1"
		"-62.7","C(C)=CC([O-])=O", // 533,"crotonate|","-CH3 : 1 | =CH- : 2 | -COO : 1 | Origin : 1 | OCCC : 1"
		"-83.39","[N+]C(CS)C([O-])=O", // 536,"cysteine|cys-L|",">CH- : 1 | -COO : 1 | -CH2- : 1 | -NH3 : 1 | -SH : 1 | Origin : 1"
		"-159.61","[N+]C(CSSCC([N+])C([O-])=O)C([O-])=O", // 537,"cystine|",">CH- : 2 | -CH2- : 2 | -COO : 2 | -NH3 : 2 | -S-S- : 1 | Origin : 1"
		"-6.19","CN1CC(=O)NC1=N", // 539,"creatinine|","ring >C= : 1 | ring -N< : 1 | =NH : 1 | ring -CH2- : 1 | -CH3 : 1 | ring >C=O : 1 | amide : 1 | ring -NH- : 1 | Origin : 1"
		"-9.55","CN", // 541,"methylamine|","methylamine : 1"
		"-1.35","[N+](C)C", // 542,"dimethylamine|","-CH3 : 2 | -NH2- : 1 | Origin : 1"
		"4.55","C[N+](C)C", // 543,"trimethylamine|","-CH3 : 3 | -NH< : 1 | Origin : 1"
		"44.52","c1ccncc1", // 544,"pyridine|","ring =CH- : 5 | ring =N- : 1 | Origin : 1 | heteroaromatic : 1"
		"-201","OC(CC([O-])=O)C([O-])=O", // 545,"malate|mal-L|",">CH- : 1 | -CH2- : 1 | -COO : 2 | -OH : 1 | Origin : 1"
		"-319.2","[H]C(=O)C(O)COP([O-])([O-])=O", // 548,"glyceraldehyde-3-phosphate|g3p|","-OPO3 : 1 | -CH2- : 1 | >CH- : 1 | -CH=O : 1 | -OH : 1 | Origin : 1"
		"-118.31","[N+]CC(=O)NCC([O-])=O", // 549,"glycylglycine|",">C=O : 1 | -CH2- : 2 | -NH3 : 1 | -COO : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-268.54","[O-]C(=O)CC(C([O-])=O)C(=O)C([O-])=O", // 551,"oxalosuccinate|",">CH- : 1 | >C=O : 1 | -CH2- : 1 | -COO : 3 | Origin : 1 | OCCO : 1"
		"-8.21","C", // 556,"ch4|methane|","methane : 1"
		"33.26","c1ccccc1", // 561,"benzene|","aromatic ring =CH- : 6 | Origin : 1 | Hydrocarbon : 1"

// proton filling not needed and not working here..  :/
//		"-1.1","[Fe+3]", // 567,"fe3|Fe3+|","Fe3 : 1"

		"-83.7","CC(C)CC([O-])=O", // 570,"isovalerate|","-CH2- : 1 | >CH- : 1 | -COO : 1 | -CH3 : 2 | Origin : 1"
		"30.6","NN", // 575,"n2h4|","N2H4 : 1"
		"-122.7","OS(S)(=O)=O", // 577,"s2o3|tsul|","S2O3H : 1"
		"-143.5","OS(=O)S(O)=O", // 578,"s2o4|","S2O4 : 1"
		"-229","OS(=O)(=O)SS(O)(=O)=O", // 579,"s3o6|","S3O6 : 1"
		"-241.82","[O-]S(=O)(=O)SSS([O-])(=O)=O", // 580,"s4o6|","-SO3 : 2 | -S-S- : 1 | Origin : 1"
		"27.62","C1CC1", // 583,"cyclopropane|","ring -CH2- : 3 | Origin : 1 | Hydrocarbon : 1 | threemember : 1"
		"206.47","CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCC1=CC(=O)c2ccccc2C1=O", // 589,"2dmmq8|2-Demethylmenaquinone-8|","aromatic ring fused to nonaromatic ring >C=  : 2 | ring >C=O : 2 | aromatic ring =CH- : 4 | ring =C< : 1 | ring =CH- : 1 | -CH2- : 15 | =CH- : 8 | =C< : 8 | -CH3 : 9 | Origin : 1 | OCCC : 4"
		"-12.6","[H]C([H])([H])Cl", // 594,"chloromethane|","chloromethane : 1"
		"-15.8","[H]C([H])(Cl)Cl", // 595,"dichloromethane|","dichloromethane : 1"
		"-15.9","[H]C(Cl)(Cl)Cl", // 596,"trichloromethane|","trichloromethane : 1"
		"-10.8","ClC(Cl)(Cl)Cl", // 597,"tetrachloromethane|","tetrachloromethane : 1"
		"-13.73","CCCl", // 598,"chloroethane|","-CH3 : 1 | -CH2- : 1 | primary -Cl : 1 | Origin : 1"
		"-18.24","ClCCCl", // 599,"1,2-dichloroethane|","-CH2- : 2 | primary -Cl : 2 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"-20.16","ClCC(Cl)Cl", // 601,"1,1,2-trichloroethane|","primary Cl2 : 2 | >CH- : 1 | primary -Cl : 1 | -CH2- : 1 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"-20.16","ClC(Cl)C(Cl)Cl", // 602,"1,1,2,2-tetrachloroethane|","primary Cl2 : 4 | >CH- : 2 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 2"
		"-13.18","CC(Cl)(Cl)Cl", // 603,"1,1,1-trichloroethane|","primary Cl3 : 3 | -CH3 : 1 | >C< : 1 | Origin : 1"
		"-17.69","ClCC(Cl)(Cl)Cl", // 604,"1,1,1,2-tetrachloroethane|","primary Cl3 : 3 | >C< : 1 | -CH2- : 1 | primary -Cl : 1 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"-13.3","ClC(Cl)(Cl)C(Cl)(Cl)Cl", // 606,"hexachloroethane|","primary Cl3 : 6 | >C< : 2 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 3"
		"-2","CCC", // 607,"propane|","-CH3 : 2 | -CH2- : 1 | Origin : 1 | Hydrocarbon : 1"
		"-0.38","CCCC", // 612,"butane|","-CH3 : 2 | -CH2- : 2 | Origin : 1 | Hydrocarbon : 1"
		"1.24","CCCCC", // 615,"pentane|","-CH2- : 3 | -CH3 : 2 | Origin : 1 | Hydrocarbon : 1"
		"8.21","ClC(Cl)=C", // 621,"1,1-dichloroethylene|","secondary Cl2 : 2 | =CH2 : 1 | =C< : 1 | Origin : 1"
		"6.52","ClC(Cl)=C(Cl)Cl", // 623,"tetrachloroethylene|","secondary Cl2 : 4 | =C< : 2 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 2"
		"24.22","Clc1ccccc1", // 624,"monochlorobenzene|","aromatic ring =CH- : 5 | aromatic ring =C< : 1 | tertiary -Cl : 1 | Origin : 1"
		"20.78","Clc1ccccc1Cl", // 625,"1,2-dichlorobenzene|","aromatic ring =CH- : 4 | aromatic ring =C< : 2 | tertiary -Cl : 2 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"18.86","Clc1cccc(Cl)c1", // 626,"1,3-dichlorobenzene|","aromatic ring =C< : 2 | tertiary -Cl : 2 | aromatic ring =CH- : 4 | Origin : 1"
		"18.86","Clc1ccc(Cl)cc1", // 627,"1,4-dichlorobenzene|","aromatic ring =C< : 2 | tertiary -Cl : 2 | aromatic ring =CH- : 4 | Origin : 1"
		"17.34","Clc1cccc(Cl)c1Cl", // 628,"1,2,3-trichlorobenzene|","aromatic ring =CH- : 3 | aromatic ring =C< : 3 | tertiary -Cl : 3 | Origin : 1 | BinaryVicinalCl : 2 | DistinctVicinalCl : 2"
		"15.42","Clc1ccc(Cl)c(Cl)c1", // 629,"1,2,4-trichlorobenzene|","aromatic ring =C< : 3 | tertiary -Cl : 3 | aromatic ring =CH- : 3 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"13.5","Clc1cc(Cl)cc(Cl)c1", // 630,"1,3,5-trichlorobenzene|","aromatic ring =C< : 3 | tertiary -Cl : 3 | aromatic ring =CH- : 3 | Origin : 1"
		"13.9","Clc1ccc(Cl)c(Cl)c1Cl", // 631,"1,2,3,4-tetrachlorobenzene|","aromatic ring =CH- : 2 | aromatic ring =C< : 4 | tertiary -Cl : 4 | Origin : 1 | BinaryVicinalCl : 3 | DistinctVicinalCl : 3"
		"11.98","Clc1cc(Cl)c(Cl)c(Cl)c1", // 632,"1,2,3,5-tetrachlorobenzene|","aromatic ring =C< : 4 | tertiary -Cl : 4 | aromatic ring =CH- : 2 | Origin : 1 | BinaryVicinalCl : 2 | DistinctVicinalCl : 2"
		"11.98","Clc1cc(Cl)c(Cl)cc1Cl", // 633,"1,2,4,5-tetrachlorobenzene|","aromatic ring =C< : 4 | tertiary -Cl : 4 | aromatic ring =CH- : 2 | Origin : 1 | BinaryVicinalCl : 2 | DistinctVicinalCl : 2"
		"10.46","Clc1cc(Cl)c(Cl)c(Cl)c1Cl", // 634,"pentachlorobenzene|","aromatic ring =CH- : 1 | aromatic ring =C< : 5 | tertiary -Cl : 5 | Origin : 1 | BinaryVicinalCl : 4 | DistinctVicinalCl : 4"
		"8.94","Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl", // 635,"hexachlorobenzene|","aromatic ring =C< : 6 | tertiary -Cl : 6 | Origin : 1 | BinaryVicinalCl : 6 | DistinctVicinalCl : 6"
		"-58.41","[O-]C(=O)c1ccccc1Cl", // 636,"2-chlorobenzoate|","aromatic ring =CH- : 4 | aromatic ring =C< : 2 | -COO : 1 | tertiary -Cl : 1 | Origin : 1 | OCCC : 1"
		"-58.41","[O-]C(=O)c1cccc(Cl)c1", // 637,"3-chlorobenzoate|","aromatic ring =CH- : 4 | aromatic ring =C< : 2 | -COO : 1 | tertiary -Cl : 1 | Origin : 1 | OCCC : 1"
		"-58.41","[O-]C(=O)c1ccc(Cl)cc1", // 638,"4-chlorobenzoate|","aromatic ring =C< : 2 | -COO : 1 | aromatic ring =CH- : 4 | tertiary -Cl : 1 | Origin : 1 | OCCC : 1"
		"-63.77","[O-]C(=O)c1ccc(Cl)cc1Cl", // 640,"2,4-dichlorobenzoate|","aromatic ring =C< : 3 | -COO : 1 | aromatic ring =CH- : 3 | tertiary -Cl : 2 | Origin : 1 | OCCC : 1"
		"-63.77","[O-]C(=O)c1c(Cl)cccc1Cl", // 642,"2,6-dichlorobenzoate|","aromatic ring =C< : 3 | -COO : 1 | aromatic ring =CH- : 3 | tertiary -Cl : 2 | Origin : 1 | OCCC : 1"
		"-61.85","[O-]C(=O)c1c([H])c([H])c(Cl)c(Cl)c1[H]", // 643,"3,4-dichlorobenzoate|","aromatic ring =CH- : 3 | aromatic ring =C< : 3 | -COO : 1 | tertiary -Cl : 2 | Origin : 1 | OCCC : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"-94.03","[O-]C(=O)c1ccccc1F", // 650,"2-fluorobenzoate|","aromatic ring =C< : 2 | -COO : 1 | aromatic ring =CH- : 4 | aromatic -F : 1 | Origin : 1 | OCCC : 1"
		"-94.03","[O-]C(=O)c1cccc(F)c1", // 651,"3-fluorobenzoate|","aromatic ring =C< : 2 | -COO : 1 | aromatic ring =CH- : 4 | aromatic -F : 1 | Origin : 1 | OCCC : 1"
		"-94.03","[O-]C(=O)c1ccc(F)cc1", // 652,"4-fluorobenzoate|","aromatic ring =C< : 2 | -COO : 1 | aromatic ring =CH- : 4 | aromatic -F : 1 | Origin : 1 | OCCC : 1"
		"-34.43","[O-]C(=O)c1ccccc1I", // 653,"2-iodobenzoate|","aromatic ring =C< : 2 | -COO : 1 | aromatic ring =CH- : 4 | aromatic -I : 1 | Origin : 1 | OCCC : 1"
		"-15.26","Oc1ccccc1Cl", // 659,"2-chlorophenol|","aromatic ring =CH- : 4 | aromatic ring =C< : 2 | -OH : 1 | tertiary -Cl : 1 | Origin : 1"
		"-15.26","Oc1cccc(Cl)c1", // 660,"3-chlorophenol|","aromatic ring =C< : 2 | -OH : 1 | aromatic ring =CH- : 4 | tertiary -Cl : 1 | Origin : 1"
		"-15.26","Oc1ccc(Cl)cc1", // 661,"4-chlorophenol|","aromatic ring =C< : 2 | -OH : 1 | aromatic ring =CH- : 4 | tertiary -Cl : 1 | Origin : 1"
		"-20.62","Oc1ccc(Cl)cc1Cl", // 663,"2,4-dichlorophenol|","aromatic ring =C< : 3 | -OH : 1 | tertiary -Cl : 2 | aromatic ring =CH- : 3 | Origin : 1"
		"-20.62","Oc1cc(Cl)ccc1Cl", // 664,"2,5-dichlorophenol|","aromatic ring =C< : 3 | -OH : 1 | aromatic ring =CH- : 3 | tertiary -Cl : 2 | Origin : 1"
		"-11.92","[O-]c1c(Cl)cccc1Cl", // 665,"2,6-dichlorophenol|","aromatic ring =C< : 3 | tertiary -Cl : 2 | aromatic ring =CH- : 3 | -Oneg : 1 | Origin : 1"
		"-24.06","Oc1cc(Cl)c(Cl)cc1Cl", // 669,"2,4,5-trichlorophenol|","aromatic ring =C< : 4 | -OH : 1 | aromatic ring =CH- : 2 | tertiary -Cl : 3 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"
		"-17.28","[O-]c1c(Cl)cc(Cl)cc1Cl", // 670,"2,4,6-trichlorophenol|","aromatic ring =C< : 4 | tertiary -Cl : 3 | aromatic ring =CH- : 2 | -Oneg : 1 | Origin : 1"
		"-18.8","[O-]c1c(Cl)c(Cl)cc(Cl)c1Cl", // 676,"2,3,5,6-tetrachlorophenol|","aromatic ring =C< : 5 | tertiary -Cl : 4 | aromatic ring =CH- : 1 | -Oneg : 1 | Origin : 1 | BinaryVicinalCl : 2 | DistinctVicinalCl : 2"
		"-20.32","[O-]c1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl", // 677,"pentachlorophenol|","aromatic ring =C< : 6 | tertiary -Cl : 5 | -Oneg : 1 | Origin : 1 | BinaryVicinalCl : 4 | DistinctVicinalCl : 4"
		"-121.26","SCCS([O-])(=O)=O", // 678,"CoM|coenzyme-M|","-SO3 : 1 | -CH2- : 2 | -SH : 1 | Origin : 1"
		"-115.39","CSCCS([O-])(=O)=O", // 679,"methyl-CoM|methyl-coenzyme-M|","-SO3 : 1 | -CH2- : 2 | -S- : 1 | -CH3 : 1 | Origin : 1"

		"-529.59","NC(=O)c1ccc[n+](c1)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC(C(O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", // 1,"NAD(+)|nad+|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 8 | ring =CH- : 6 | ring =N- : 3 | ring -O- : 2 | ring =C< : 2 | -OH : 4 | -NH2 : 2 | -CH2- : 2 | -OPO3- : 1 | -OPO2- : 1 | dbl sgl ring =N< : 1 | >C=O : 1 | amide : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 3"
		"-524.32","NC(=O)C1=CN(C=CC1)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC(C(O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", // 2,"NADH|nadh|","ring -N< : 2 | two fused rings >C= : 2 | ring -CH< : 8 | ring =CH- : 5 | ring =N- : 3 | ring -O- : 2 | ring =C< : 2 | -OH : 4 | -NH2 : 2 | -CH2- : 2 | -OPO3- : 1 | -OPO2- : 1 | ring -CH2- : 1 | >C=O : 1 | amide : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 2 | NADH : 1"
		"-736.82","NC(=O)C1=CN(C=CC1)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC(C(OP([O-])([O-])=O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", // 3,"NADPH|","ring -CH< : 8 | ring -N< : 2 | ring -O- : 2 | two fused rings >C= : 2 | ring =CH- : 5 | -OPO3 : 1 | ring =N- : 3 | -OH : 3 | -CH2- : 2 | ring =C< : 2 | -OPO3- : 1 | -NH2 : 2 | -OPO2- : 1 | ring -CH2- : 1 | >C=O : 1 | amide : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 2 | NADH : 1"
		"-742.09","NC(=O)c1ccc[n+](c1)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC(C(OP([O-])([O-])=O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", // 4,"NADP(+)|","ring -CH< : 8 | ring -N< : 1 | ring -O- : 2 | two fused rings >C= : 2 | ring =CH- : 6 | -OPO3 : 1 | ring =N- : 3 | -OH : 3 | -CH2- : 2 | ring =C< : 2 | -OPO3- : 1 | -NH2 : 2 | -OPO2- : 1 | dbl sgl ring =N< : 1 | >C=O : 1 | amide : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 3"
		"-751.99","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n2cnc3c(N)ncnc23)C(O)C(=O)NCCC(=O)NCCS", // 7,"Coenzyme A|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 6 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 2 | -SH : 1 | amide : 2 | -NH- : 2 | Origin : 1 | heteroaromatic : 2 | COA : 1"
		"-605.77","C(OP(=O)([O-])OP(=O)([O-])[O-])C1(OC(C(O)C1(O))N2(C=CC(=O)NC2(=O)))", // 10,"UDP|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 2 | ring =CH- : 2 | -OH : 2 | -CH2- : 1 | -OPO2- : 1 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"

// only 2 heteroaromatic rings by jankowski
//		"-538.18","Cc1cc2Nc3c([nH]c(=O)[nH]c3=O)N(CC(O)C(O)C(O)COP([O-])(=O)OP([O-])(=O)OCC3OC(C(O)C3O)n3cnc4c(N)ncnc34)c2cc1C", // 11,"FADH2|","ring -N< : 2 | two fused rings >C= : 4 | aromatic ring fused to nonaromatic ring >C=  : 2 | -CH2- : 3 | aromatic ring =CH- : 2 | >CH- : 3 | ring -NH- : 3 | ring >C=O : 2 | aromatic ring =C< : 2 | -OH : 5 | -CH3 : 2 | -OPO3- : 1 | -OPO2- : 1 | ring -CH< : 4 | ring -O- : 1 | ring =CH- : 2 | ring =N- : 3 | ring =C< : 1 | -NH2 : 1 | amide : 2 | Origin : 1 | OCCC : 1 | heteroaromatic : 2"

		"-257.85","Nc1ncnc2n(cnc12)C1OC(COP([O-])([O-])=O)C(O)C1O", // 13,"AMP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring =C< : 1 | -OH : 2 | -NH2 : 1 | -CH2- : 1 | -OPO3 : 1 | Origin : 1 | heteroaromatic : 2"
		"-76.11","[N+]C(CCSCC1OC(C(O)C1O)n1cnc2c(N)ncnc12)C([O-])=O", // 14,"S-adenosyl-L-homocysteine|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring =C< : 1 | -OH : 2 | -NH2 : 1 | -CH2- : 3 | -S- : 1 | >CH- : 1 | -COO : 1 | -NH3 : 1 | Origin : 1 | heteroaromatic : 2"
		"-785.83","CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 16,"acetyl-CoA|accoa|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 6 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 3 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"

// only 1 heteroaromatic ring by jankowski
// no OCCC match -> double bond in larger heteroring?
		"-516.3","N1=C(N)NC(=O)c2ncn(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)c21", // 19,"GDP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 1 | ring =N- : 2 | ring -O- : 1 | ring >C=O : 1 | ring =C< : 1 | -OH : 2 | -NH2 : 1 | -CH2- : 1 | -OPO2- : 1 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-756.9","C(C1(C(C(C(O1)N2(C=CC(NC2=O)=O))O)O))OP(OP(OC3(C(C(C(C(O3)CO)O)O)NC(C)=O))([O-])=O)([O-])=O", // 24,"UDP-N-acetyl-D-glucosamine|","ring -CH< : 9 | ring -O- : 2 | ring >C=O : 2 | ring =CH- : 2 | -OH : 5 | -CH2- : 2 | -OPO3- : 1 | -OPO2- : 1 | >C=O : 1 | -CH3 : 1 | amide : 3 | ring -N< : 1 | ring -NH- : 1 | -NH- : 1 | Origin : 1 | OCCC : 1"

// only 1 heteroaromatic ring by jankowski
//		"-724.3","Nc1nc2n(cnc2c(=O)[nH]1)C1OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C1O", // 25,"GTP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 1 | ring =N- : 2 | ring -O- : 1 | ring >C=O : 1 | ring =C< : 1 | -OH : 2 | -NH2 : 1 | -CH2- : 1 | -OPO2- : 2 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-165.75","[NH3+]C(CC([O-])=O)C([O-])=O", // 27,"L-aspartate|asp-L|",">CH- : 1 | -CH2- : 1 | -COO : 2 | -NH3 : 1 | Origin : 1"
		"-220.95","[N+]C(CCC(=O)NC(CS)C(=O)NCC([O-])=O)C([O-])=O", // 28,"reduced-glutathione|",">CH- : 2 | >C=O : 2 | -CH2- : 4 | -SH : 1 | -COO : 2 | -NH3 : 1 | amide : 2 | -NH- : 2 | Origin : 1"
		"-762.55","C2(C(O)C(O)C(COP(=O)([O-])OP(=O)(OC1(OC(CO)C(O)C(O)C1(O)))[O-])O2)N3(C=CC(=O)NC3(=O))", // 29,"UDPgalactose|UDPglucose|","ring -CH< : 9 | ring -O- : 2 | ring >C=O : 2 | ring =CH- : 2 | -OH : 6 | -CH2- : 2 | -OPO3- : 1 | -OPO2- : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-345.38","C(C2(C(C(C(N1(C(N=C(C=C1)N)=O))O2)O)O))OP([O-])([O-])=O", // 31,"CMP|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -CH2- : 1 | ring =C< : 1 | -OPO3 : 1 | -NH2 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-346.33","C12(N(CC(O)C(O)C(O)COP([O-])(=O)[O-])c3(c(NC=2(C(=O)NC(=O)N1))cc(C)c(C)c3))", // 34,"FMNH2|fmnred|","ring -N< : 1 | aromatic ring fused to nonaromatic ring >C=  : 2 | two fused rings >C= : 2 | -CH2- : 2 | aromatic ring =CH- : 2 | >CH- : 3 | ring -NH- : 3 | aromatic ring =C< : 2 | ring >C=O : 2 | -OH : 3 | -CH3 : 2 | -OPO3 : 1 | amide : 2 | Origin : 1 | OCCC : 1"
		"-73.77","[N+]C(CCCNC(N)=[N+])C([O-])=O", // 35,"L-arginine|arg|","-CH2- : 3 | >CH- : 1 | -COO : 1 | -NH3 : 1 | -NH- : 1 | =C< : 1 | -NH2 : 1 | =NH2 : 1 | Origin : 1"
		"-761.38","C(OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])C1(OC(C(O)C1(O))N2(C=CC(N)=NC2(=O)))", // 36,"CTP|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -CH2- : 1 | ring =C< : 1 | -OPO2- : 2 | -NH2 : 1 | -OPO3 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-121.69","[N+]C(CCC(N)=O)C([O-])=O", // 37,"L-glutamine|gln-L|",">CH- : 1 | -CH2- : 2 | -COO : 1 | -NH3 : 1 | >C=O : 1 | amide : 1 | -NH2 : 1 | Origin : 1"
		"-75.91","CSCCC([N+])C([O-])=O", // 39,"L-methionine|met-L|",">CH- : 1 | -COO : 1 | -CH2- : 2 | -NH3 : 1 | -S- : 1 | -CH3 : 1 | Origin : 1"
		"-813.77","C(OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])C1(OC(C(O)C1(O))N2(C=CC(=O)NC2(=O)))", // 41,"UTP|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 2 | ring =CH- : 2 | -OH : 2 | -CH2- : 1 | -OPO2- : 2 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-25","[N+]C(Cc1c[nH]c2ccccc12)C([O-])=O", // 42,"L-tryptophan|trp-L|","aromatic ring fused to nonaromatic ring >C=  : 2 | ring =C< : 1 | aromatic ring =CH- : 4 | -CH2- : 1 | ring =CH- : 1 | ring -NH- : 1 | >CH- : 1 | -COO : 1 | -NH3 : 1 | Origin : 1 | heteroaromatic : 1"
		"-51.05","[N+]C(Cc1ccccc1)C([O-])=O", // 43,"L-phenylalanine|phe-L|","aromatic ring =C< : 1 | -CH2- : 1 | aromatic ring =CH- : 5 | >CH- : 1 | -COO : 1 | -NH3 : 1 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"-729.58","OC1C(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)OC(C1O)n1cnc2c1nc[nH]c2=O", // 45,"ITP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 2 | ring -O- : 1 | ring >C=O : 1 | -OH : 2 | -CH2- : 1 | -OPO2- : 2 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-90.53","[N+]C(Cc1ccc(O)cc1)C([O-])=O", // 46,"L-tyrosine|tyr-L|",">CH- : 1 | -COO : 1 | -CH2- : 1 | -NH3 : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -OH : 1 | Origin : 1"
		"-863.66","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCSC(=O)CC([O-])=O", // 47,"malonyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 7 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 3 | -COO : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-431.67","C(OP(=O)([O-])[O-])C1(OC(O)(CO)C(O)C1(O))", // 49,"D-fructose 6-phosphate|","ring -CH< : 3 | ring -O- : 1 | -CH2- : 2 | ring >C< : 1 | -OH : 4 | -OPO3 : 1 | Origin : 1"
		"-377.65","OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O", // 51,"sucrose|sucr|","ring >C< : 1 | -O- : 1 | ring -O- : 2 | ring -CH< : 8 | -CH2- : 3 | -OH : 8 | Origin : 1"
		"-862.04","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCSC(=O)CCC([O-])=O", // 52,"succinyl-CoA|succoa|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 8 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 3 | -COO : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"

// only 1 heteroaromatic ring by jankowski
//		"-521.58","OC1C(COP([O-])(=O)OP([O-])([O-])=O)OC(C1O)n1cnc2c1nc[nH]c2=O", // 56,"IDP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 2 | ring -O- : 1 | ring >C=O : 1 | -OH : 2 | -CH2- : 1 | -OPO2- : 1 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-397.77","C(OP([O-])([O-])=O)C1OC(N2C=CC(=O)NC2(=O))C(O)C1(O)", // 57,"UMP|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 2 | ring =CH- : 2 | -OH : 2 | -CH2- : 1 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-553.38","C(C2(C(C(C(N1(C(N=C(C=C1)N)=O))O2)O)O))OP(OP([O-])([O-])=O)([O-])=O", // 60,"CDP|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -CH2- : 1 | ring =C< : 1 | -OPO2- : 1 | -NH2 : 1 | -OPO3 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-394.12","C(OP([O-])(=O)[O-])C1(C(O)C(O)C(O)O1)", // 63,"D-ribose 5-phosphate|","-OPO3 : 1 | -CH2- : 1 | ring -CH< : 4 | ring -O- : 1 | -OH : 3 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"-313.58","OC1C(COP([O-])([O-])=O)OC(C1O)n1cnc2c1nc[nH]c2=O", // 69,"IMP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 2 | ring -O- : 1 | ring >C=O : 1 | -OH : 2 | -CH2- : 1 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-634.01","Nc1ncnc2n(cnc12)C1CC(O)C(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)O1", // 70,"dATP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 3 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring -CH2- : 1 | ring =C< : 1 | -NH2 : 1 | -CH2- : 1 | -OH : 1 | -OPO2- : 2 | -OPO3 : 1 | Origin : 1 | heteroaromatic : 2"
		"-45.63","[N+]C(Cc1c[nH]cn1)C([O-])=O", // 73,"L-histidine|","ring =C< : 1 | -CH2- : 1 | ring =CH- : 2 | ring =N- : 1 | >CH- : 1 | ring -NH- : 1 | -COO : 1 | -NH3 : 1 | Origin : 1 | heteroaromatic : 1"
		"-219.96","OC1C(O)C(O)C(O)C(O)C1O", // 74,"myo-inositol|","ring -CH< : 6 | -OH : 6 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"-308.3","Nc1nc2n(cnc2c(=O)[nH]1)C1OC(COP([O-])([O-])=O)C(O)C1O", // 75,"GMP|","two fused rings >C= : 2 | ring -N< : 1 | ring =N- : 2 | ring -CH< : 4 | ring =CH- : 1 | ring >C=O : 1 | ring =C< : 1 | ring -O- : 1 | -NH2 : 1 | -OH : 2 | -CH2- : 1 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-62.54","[O-]C(=O)C1CCCN1", // 78,"S-proline|","ring -CH< : 1 | -COO : 1 | ring -CH2- : 3 | ring -NH- : 1 | Origin : 1"
		"-123.31","[N+]C(CC(N)=O)C([O-])=O", // 79,"L-asparagine|asn-L|",">CH- : 1 | -CH2- : 1 | -COO : 1 | -NH3 : 1 | >C=O : 1 | amide : 1 | -NH2 : 1 | Origin : 1"
		"-763.15","CCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 81,"palmitoyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 20 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 3 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-81.77","[N+]C(CCS)C([O-])=O", // 82,"L-homocysteine|",">CH- : 1 | -CH2- : 2 | -COO : 1 | -NH3 : 1 | -SH : 1 | Origin : 1"
		"-86.49","CC(C)C([N+])C([O-])=O", // 89,"L-valine|val-L|",">CH- : 2 | -COO : 1 | -NH3 : 1 | -CH3 : 2 | Origin : 1"
		"-124.34","CC(O)C([N+])C([O-])=O", // 91,"L-threonine|thr-L|",">CH- : 2 | -COO : 1 | -NH3 : 1 | -CH3 : 1 | -OH : 1 | Origin : 1"
		"-216.02","OCC1OC(=O)C(O)C(O)C1O", // 93,"D-glucono-1,5-lactone|","ring -CH< : 4 | ring -O-CO- : 1 | -CH2- : 1 | -OH : 4 | Origin : 1"
		"-426.01","Nc1ncnc2n(cnc12)C1CC(O)C(COP([O-])(=O)OP([O-])([O-])=O)O1", // 94,"dADP|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 3 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring -CH2- : 1 | ring =C< : 1 | -NH2 : 1 | -CH2- : 1 | -OH : 1 | -OPO2- : 1 | -OPO3 : 1 | Origin : 1 | heteroaromatic : 2"
		"-45.35","Nc1ncnc2n(cnc12)C1OC(CO)C(O)C1O", // 96,"adenosine|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring =C< : 1 | -OH : 3 | -NH2 : 1 | -CH2- : 1 | Origin : 1 | heteroaromatic : 2"
		"-164.13","[NH3+]C(CCC([O-])=O)C([O-])=O", // 97,"D-glutamate|glu-L|",">CH- : 1 | -CH2- : 2 | -COO : 2 | -NH3 : 1 | Origin : 1"

// wrong SMILES -> no -COO group by jankowski
//		"-586.8","OC(=O)C(COP([O-])([O-])=O)OP([O-])([O-])=O", // 103,"2,3-bisphospho-D-glyceric acid|1,3-bpg|","-CO-OPO3- : 1 | >CH- : 1 | -CH2- : 1 | -OH : 1 | -OPO3 : 1 | Origin : 1"

		"-305.54","C(C2(C(CC(N1(C(N=C(C=C1)N)=O))O2)O))OP([O-])([O-])=O", // 104,"dCMP|","ring -CH< : 3 | ring >C=O : 1 | ring =CH- : 2 | ring -O- : 1 | ring -CH2- : 1 | ring =N- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 1 | -NH2 : 1 | -OPO3 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-169.49","OC1C=CC(=CC1OC(=C)C([O-])=O)C([O-])=O", // 108,"chorismate|","ring -CH< : 2 | ring =CH- : 3 | -O- : 1 | ring =C< : 1 | -OH : 1 | =C< : 1 | -COO : 2 | =CH2 : 1 | Origin : 1 | OCCC : 2 | CCCC : 1"
		"-268.66","OCC(O)C(O)C(O)C(O)C([O-])=O", // 110,"D-gluconate|glcn|",">CH- : 4 | -OH : 5 | -COO : 1 | -CH2- : 1 | Origin : 1"
		"-122.53","[N+]C(CCO)C([O-])=O", // 114,"L-homoserine|",">CH- : 1 | -CH2- : 2 | -COO : 1 | -NH3 : 1 | -OH : 1 | Origin : 1"
		"-326.4","[H]C1(OC(O)(CC(O)C1NC(C)=O)C([O-])=O)C(O)C(O)CO", // 116,"N-acetylneuraminate|","ring -CH< : 3 | ring -O- : 1 | >CH- : 2 | ring >C< : 1 | -OH : 5 | ring -CH2- : 1 | >C=O : 1 | -COO : 1 | -CH2- : 1 | -CH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-209.4","OCC(=O)C(O)C(O)C(=O)CO", // 117,"5-dehydro-D-fructose|",">CH- : 2 | >C=O : 2 | -OH : 4 | -CH2- : 2 | Origin : 1"
		"-355.62","[H]C(=O)C(O)C(O)COP([O-])([O-])=O", // 118,"D-erythrose-4-phosphate|","-OPO3 : 1 | -CH2- : 1 | >CH- : 2 | -OH : 2 | -CH=O : 1 | Origin : 1"
		"-18.43","[H]C12CCC3=CC(=O)CCC3(C)C1([H])CCC1(C)C(=O)CCC21[H]", // 119,"4-androstene-3,17-dione|","two fused rings >CH- : 3 | ring -CH2- : 8 | two fused rings >C< : 2 | two fused rings >C= : 1 | -CH3 : 2 | ring >C=O : 2 | ring =CH- : 1 | Origin : 1 | OCCC : 1"

// only 1 heteroaromatic ring by jankowski
//		"-101.08","OCC1OC(C(O)C1O)n1cnc2c(O)ncnc12", // 121,"inosine|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 2 | ring -O- : 1 | ring >C=O : 1 | -OH : 3 | -CH2- : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-185.27","OCC1OC(C(O)C1O)N2C=CC(=O)NC2=O", // 124,"uridine|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 2 | ring =CH- : 2 | -OH : 3 | -CH2- : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-479.69","C(OP([O-])(=O)OP(OCC[N+](C)(C)C)([O-])=O)C1(OC(C(O)C1(O))N2(C=CC(=NC2=O)N))", // 126,"CDP-choline|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -CH2- : 3 | ring =C< : 1 | -OPO3- : 1 | -NH2 : 1 | -OPO2- : 1 | >N< : 1 | -CH3 : 3 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-862.82","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCSC(=O)C([O-])=O", // 128,"oxalyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 6 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 3 | -COO : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | OCCO : 1 | heteroaromatic : 2"
		"-123.87","[N+]C(CCCNC(N)=O)C([O-])=O", // 129,"L-citrulline|citr-L|","-CH2- : 3 | >CH- : 1 | -COO : 1 | -NH3 : 1 | >C=O : 1 | amide : 2 | -NH- : 1 | -NH2 : 1 | Origin : 1"
		"-812.61","CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 131,"acetoacetyl-CoA|aacoa|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 7 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 4 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-261.5","OC1OC(C(O)C(O)C1O)C([O-])=O", // 132,"α-D-galacturonate|CHEBI:58658","ring -CH< : 5 | ring -O- : 1 | -COO : 1 | -OH : 4 | Origin : 1"
		"-151.52","[O-]C(=O)C1CC(=O)NC(=O)N1", // 133,"S-dihydroorotate|","ring -CH< : 1 | ring -CH2- : 1 | -COO : 1 | ring >C=O : 2 | amide : 2 | ring -NH- : 2 | Origin : 1"
		"-481.16","OC(COP([O-])([O-])=O)C(O)C(O)C(O)C([O-])=O", // 136,"6-phospho-D-gluconate|","-OPO3 : 1 | -CH2- : 1 | >CH- : 4 | -OH : 4 | -COO : 1 | Origin : 1"
		"-395.53","C(OP(=O)([O-])[O-])C1(OC(O)C([N+])C(O)C1(O))", // 138,"D-glucosamine-6-phosphate|","-OPO3 : 1 | -CH2- : 1 | ring -CH< : 5 | ring -O- : 1 | -OH : 3 | -NH3 : 1 | Origin : 1"
		"-644.17","C(C1(C(C(C(O1)(COP([O-])([O-])=O)O)O)O))OP([O-])([O-])=O", // 139,"D-fructose-1,6-bisphosphate|","ring >C< : 1 | ring -O- : 1 | ring -CH< : 3 | -CH2- : 2 | -OH : 3 | -OPO3 : 2 | Origin : 1"
		"-376.76","OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O", // 140,"laminaribiose|beta-laminaribiose|","ring -CH< : 10 | -O- : 1 | -OH : 8 | ring -O- : 2 | -CH2- : 2 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"-476.46","Nc1nc2n(cnc2c(=O)[nH]1)C1CC(O)C(COP([O-])(=O)OP([O-])([O-])=O)O1", // 141,"dGDP|","two fused rings >C= : 2 | ring -N< : 1 | ring =N- : 2 | ring -CH< : 3 | ring =CH- : 1 | ring >C=O : 1 | ring =C< : 1 | ring -O- : 1 | ring -CH2- : 1 | -NH2 : 1 | -CH2- : 1 | -OH : 1 | -OPO2- : 1 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"
//		"-268.46","Nc1nc2n(cnc2c(=O)[nH]1)C1CC(O)C(COP([O-])([O-])=O)O1", // 142,"dGMP|","two fused rings >C= : 2 | ring -N< : 1 | ring =N- : 2 | ring -CH< : 3 | ring =CH- : 1 | ring >C=O : 1 | ring =C< : 1 | ring -O- : 1 | ring -CH2- : 1 | -NH2 : 1 | -CH2- : 1 | -OH : 1 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-566.34","CC1(=CN(C(=O)NC1(=O))C2(CC(O)C(COP(=O)([O-])OP(=O)([O-])[O-])O2))", // 143,"dTDP|","ring -CH< : 3 | ring -O- : 1 | ring -CH2- : 1 | ring >C=O : 2 | ring =CH- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 1 | -CH3 : 1 | -OPO2- : 1 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-358.34","CC1(=CN(C(=O)NC1(=O))C2(CC(O)C(COP(=O)([O-])[O-])O2))", // 144,"dTMP|","ring -CH< : 3 | ring -O- : 1 | ring -CH2- : 1 | ring >C=O : 2 | ring =CH- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 1 | -CH3 : 1 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"

// only 1 heteroaromatic ring by jankowski
//		"-95.8","Nc1nc2n(cnc2c(=O)[nH]1)C1OC(CO)C(O)C1O", // 149,"guanosine|","two fused rings >C= : 2 | ring -N< : 1 | ring =N- : 2 | ring -CH< : 4 | ring =CH- : 1 | ring >C=O : 1 | ring =C< : 1 | ring -O- : 1 | -NH2 : 1 | -OH : 3 | -CH2- : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-115.26","CC(=O)NC(CCC[N+])C([O-])=O", // 153,"N2-acetyl-L-ornithine|CHEBI:16543",">CH- : 1 | -COO : 1 | -CH2- : 3 | >C=O : 1 | -CH3 : 1 | -NH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-203.96","NC(=O)NC(CC([O-])=O)C([O-])=O", // 154,"N-carbamoyl-L-aspartate|",">CH- : 1 | -CH2- : 1 | -COO : 2 | >C=O : 1 | amide : 2 | -NH- : 1 | -NH2 : 1 | Origin : 1"
		"-157.18","[H]C(=[N+])NC(CCC(O)=O)C([O-])=O", // 155,"N-formimino-L-glutamate|",">CH- : 1 | -CH2- : 2 | -NH- : 1 | -COO : 2 | =CH- : 1 | =NH2 : 1 | Origin : 1"
		"-113.05","[H]C(=O)CC([N+])C([O-])=O", // 156,"L-aspartate-4-semialdehyde|",">CH- : 1 | -CH2- : 1 | -COO : 1 | -NH3 : 1 | -CH=O : 1 | Origin : 1"
		"-337.74","NC(=O)c1ccc[n+](c1)C1OC(COP(O)([O-])=O)C(O)C1O", // 157,"beta-nicotinamide-mononucleotide|CHEBI:16171","ring -CH< : 4 | dbl sgl ring =N< : 1 | ring -O- : 1 | ring =CH- : 4 | -OH : 2 | ring =C< : 1 | -CH2- : 1 | >C=O : 1 | -OPO3 : 1 | amide : 1 | -NH2 : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 1"
		"-774.34","CC1(=CN(C(=O)NC1(=O))C2(CC(O)C(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O2))", // 158,"dTTP|","ring -CH< : 3 | ring -O- : 1 | ring -CH2- : 1 | ring >C=O : 2 | ring =CH- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 1 | -CH3 : 1 | -OPO2- : 2 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-16.98","[H]C12CCC3(C)C(=O)CCC3([H])C1([H])CCc1cc(O)ccc21", // 159,"estrone|","two fused rings >CH- : 3 | two fused rings >C< : 1 | ring -CH2- : 6 | ring >C=O : 1 | -CH3 : 1 | aromatic ring fused to nonaromatic ring >C=  : 2 | aromatic ring =CH- : 3 | aromatic ring =C< : 1 | -OH : 1 | Origin : 1"
		"-132.88","C1(N=C(C=CN1C2(C(O)C(O)C(CO)O2))N)(=O)", // 161,"cytidine|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 3 | ring =N- : 1 | -CH2- : 1 | ring =C< : 1 | -NH2 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-536.13","OCC1OC(OCC2OC(OC3(CO)OC(CO)C(O)C3O)C(O)C(O)C2O)C(O)C(O)C1O", // 165,"raffinose|","ring >C< : 1 | -O- : 2 | ring -O- : 3 | ring -CH< : 13 | -CH2- : 4 | -OH : 11 | Origin : 1"
		"-171.29","OC1CC(=CC(O)C1O)C([O-])=O", // 166,"shikimate|","ring =C< : 1 | ring -CH2- : 1 | ring =CH- : 1 | -COO : 1 | ring -CH< : 3 | -OH : 3 | Origin : 1 | OCCC : 1"
		"-622.63","Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OC2OC(CO)C(O)C(O)C2O)C(O)C1O", // 167,"ADPglucose|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 9 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 2 | ring =C< : 1 | -OH : 6 | -NH2 : 1 | -CH2- : 2 | -OPO3- : 1 | -OPO2- : 1 | Origin : 1 | heteroaromatic : 2"
		"-710.16","C(C2(C(C(C(N1(C(N=C(C=C1)N)=O))O2)O)O))OP(OP(OC3(C(C(C(C(O3)CO)O)O)O))([O-])=O)([O-])=O", // 169,"CDPglucose|CHEBI:58660","ring -CH< : 9 | ring -O- : 2 | ring >C=O : 1 | ring =CH- : 2 | -OH : 6 | ring =N- : 1 | -CH2- : 2 | ring =C< : 1 | -OPO3- : 1 | -NH2 : 1 | -OPO2- : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-189.02","OCC(O)C(O)C(O)CO", // 172,"arabitol|xylitol|ribitol|",">CH- : 3 | -OH : 5 | -CH2- : 2 | Origin : 1"
		"-83.86","NC(C(O)=O)c1ccc(O)cc1", // 174,"D-4-hydroxyphenylglycine|","-COO : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | >CH- : 1 | -NH2 : 1 | -OH : 1 | Origin : 1"
		"-239.04","OC(C(O)C([O-])=O)C([O-])=O", // 175,"meso-tartrate|",">CH- : 2 | -COO : 2 | -OH : 2 | Origin : 1"
		"-536.39","C(C2(C(C(C(N1(C(N=C(C=C1)N)=O))O2)O)O))OP(OP(OCC[N+])([O-])=O)([O-])=O", // 178,"CDP-ethanolamine|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -CH2- : 3 | ring =C< : 1 | -OPO3- : 1 | -NH2 : 1 | -OPO2- : 1 | -NH3 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-191.87","Nc1ncnc2n(cnc12)C1OC2COP([O-])(=O)OC2C1O", // 179,"adenosine-3'5'-cyclicphosphate|CHEBI:17489","ring -N< : 1 | ring -CH< : 2 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 2 | ring =N- : 3 | two fused rings >CH- : 2 | -OH : 1 | ring =C< : 1 | ring -OPO2- : 1 | ring -CH2- : 1 | -NH2 : 1 | Origin : 1 | heteroaromatic : 2"
		"-135.27","CCCCCC(O)C=CC1C(O)CC(=O)C1CCCCCCC([O-])=O", // 184,"prostaglandin E1|","ring -CH< : 3 | ring >C=O : 1 | -CH2- : 10 | =CH- : 2 | ring -CH2- : 1 | -OH : 2 | >CH- : 1 | -COO : 1 | -CH3 : 1 | Origin : 1"
		"-112.91","CCCCCC(O)C=CC1C(O)CC(=O)C1CC=CCCCC([O-])=O", // 185,"prostaglandin E2|","ring -CH< : 3 | ring >C=O : 1 | -CH2- : 8 | =CH- : 4 | ring -CH2- : 1 | -OH : 2 | >CH- : 1 | -COO : 1 | -CH3 : 1 | Origin : 1"
		"-371.9","OC(COP([O-])([O-])=O)C([O-])=O", // 188,"3-phospho-D-glycerate|pg3|CHEBI:17794","-OPO3 : 1 | -CH2- : 1 | >CH- : 1 | -COO : 1 | -OH : 1 | Origin : 1"
		"-163.98","NC(=O)NC(O)C([O-])=O", // 189,"--ureidoglycolate|CHEBI:57296",">CH- : 1 | -COO : 1 | -OH : 1 | >C=O : 1 | amide : 2 | -NH- : 1 | -NH2 : 1 | Origin : 1"
		"-394.12","OCC1OC(OP([O-])([O-])=O)C(O)C1O", // 190,"D-ribose-1-phosphate|","ring -CH< : 4 | ring -O- : 1 | -OPO3 : 1 | -OH : 3 | -CH2- : 1 | Origin : 1"
		"-193.73","CC(=O)NC(CCC([O-])=O)C([O-])=O", // 191,"N-acetyl-L-glutamate|",">CH- : 1 | -CH2- : 2 | -COO : 2 | >C=O : 1 | -CH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-371.9","[H]C(CO)(OP([O-])([O-])=O)C([O-])=O", // 192,"2-phospho-D-glycerate|pg2|CHEBI:17835","-OPO3 : 1 | >CH- : 1 | -COO : 1 | -CH2- : 1 | -OH : 1 | Origin : 1"
		"-212.63","C(O)C1(OC(O)C(NC(C)=O)C(O)C1(O))", // 194,"N-acetyl-D-mannosamine|N-acetyl-D-glucosamine|","ring -CH< : 5 | -OH : 4 | ring -O- : 1 | >C=O : 1 | -CH3 : 1 | -CH2- : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-643.28","C(C1(OC(C(C(C1O)O)O)OP(=O)([O-])[O-]))OP(=O)([O-])[O-]", // 195,"D-glucose 1,6-bisphosphate|","-OPO3 : 2 | ring -CH< : 5 | ring -O- : 1 | -OH : 3 | -CH2- : 1 | Origin : 1"
		"-354.28","OCC1OC(CC1O)OP([O-])([O-])=O", // 196,"2-deoxy-alpha-D-ribose-1-phosphate|CHEBI:11563","ring -CH< : 3 | ring -O- : 1 | ring -CH2- : 1 | -OPO3 : 1 | -CH2- : 1 | -OH : 2 | Origin : 1"
		"-354.28","C(OP([O-])(=O)[O-])C1(OC(O)CC1(O))", // 197,"2-deoxy-D-ribose-5-phosphate|CHEBI:16132","-OPO3 : 1 | -CH2- : 1 | ring -CH< : 3 | ring -O- : 1 | ring -CH2- : 1 | -OH : 2 | Origin : 1"
		"-680.33","CC1(=CN(C(=O)NC1(=O))C3(CC(O)C(COP(=O)([O-])OP(=O)([O-])OC2(OC(C)C(=O)C(O)C2(O)))O3))", // 198,"dTDP-4-dehydro-6-deoxy-D-glucose|CHEBI:16128","ring -CH< : 7 | ring -O- : 2 | ring -CH2- : 1 | ring >C=O : 3 | ring =CH- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 3 | -CH3 : 2 | -OPO3- : 1 | -OPO2- : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-513.54","C(C2(C(CC(N1(C(N=C(C=C1)N)=O))O2)O))OP(OP(=O)([O-])[O-])([O-])=O", // 199,"dCDP|","ring -CH< : 3 | ring >C=O : 1 | ring =CH- : 2 | ring -O- : 1 | ring -CH2- : 1 | ring =N- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 1 | -NH2 : 1 | -OPO2- : 1 | -OPO3 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-784.18","[H]C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 203,"formyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 6 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 2 | -CH=O : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-723.12","CC1(=CN(C(=O)NC1(=O))C3(CC(O)C(COP(=O)([O-])OP(=O)([O-])OC2(OC(CO)C(O)C(O)C2(O)))O3))", // 207,"dTDPglucose|","ring -CH< : 8 | ring -O- : 2 | ring -CH2- : 1 | ring >C=O : 2 | ring =CH- : 1 | ring =C< : 1 | -CH2- : 2 | -OH : 5 | -CH3 : 1 | -OPO3- : 1 | -OPO2- : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"
		"-761.78","CC=CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 209,"cis-but-2-enoyl-CoA|trans-but-2-enoyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 6 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 3 | =CH- : 2 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 2"
		"-539.49","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS", // 210,"3'-dephospho-CoA|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring =C< : 1 | -OH : 3 | -NH2 : 1 | -CH2- : 6 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 2 | -SH : 1 | amide : 2 | -NH- : 2 | Origin : 1 | heteroaromatic : 2"
		"-169.49","OC1C(OC(=C)C([O-])=O)C=CC=C1C([O-])=O", // 211,"isochorismate|","ring -CH< : 2 | -O- : 1 | ring =CH- : 3 | ring =C< : 1 | -OH : 1 | =C< : 1 | -COO : 2 | =CH2 : 1 | Origin : 1 | OCCC : 2 | CCCC : 1"
		"-965.77","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCSSCC(NC(=O)CCC([N+])C(O)=O)C(=O)NCC([O-])=O", // 215,"CoA-glutathione|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 10 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 3 | -CH3 : 2 | >C=O : 4 | -S-S- : 1 | -COO : 2 | -NH3 : 1 | amide : 4 | -NH- : 4 | Origin : 1 | heteroaromatic : 2"

// no heteroaromatic ring by jankowski
//		"-280.96","Nc1ccn(C2OC3COP([O-])(=O)OC3C2O)c(=O)n1", // 218,"cytidine-2'3'-cyclicphosphate|CHEBI:17065","ring -CH< : 2 | two fused rings >CH- : 2 | ring -O- : 2 | ring >C=O : 1 | ring =CH- : 2 | ring =N- : 1 | ring -OPO2- : 1 | -CH2- : 1 | ring =C< : 1 | -OH : 1 | -NH2 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"

		"-23.54","[H]C12CCC3(C)C(O)CCC3([H])C1([H])CCc1cc(O)ccc21", // 219,"estradiol-17beta|CHEBI:16469","two fused rings >CH- : 3 | two fused rings >C< : 1 | ring -CH2- : 6 | ring -CH< : 1 | -CH3 : 1 | aromatic ring fused to nonaromatic ring >C=  : 2 | -OH : 2 | aromatic ring =CH- : 3 | aromatic ring =C< : 1 | Origin : 1"
		"-162.51","[NH3+]C(CCCC([O-])=O)C([O-])=O", // 220,"L-2-aminoadipate|","-CH2- : 3 | >CH- : 1 | -COO : 2 | -NH3 : 1 | Origin : 1"
		"-161.6","CC(=O)OCC([N+])C([O-])=O", // 222,"acetyl-L-serine|CHEBI:17981",">CH- : 1 | -CH2- : 1 | -COO : 1 | -NH3 : 1 | -O-CO- : 1 | -CH3 : 1 | Origin : 1"
		"-376.76","OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O", // 227,"alpha,alpha-trehalose|","ring -CH< : 10 | -O- : 1 | ring -O- : 2 | -OH : 8 | -CH2- : 2 | Origin : 1"
		"-431.67","C(O)C1(OC(O)(COP(=O)([O-])[O-])C(O)C1(O))", // 229,"fructose-1-phosphate|CHEBI:18105","ring >C< : 1 | ring -CH< : 3 | ring -O- : 1 | -CH2- : 2 | -OH : 4 | -OPO3 : 1 | Origin : 1"
		"-437.94","OCC(O)C(O)C(O)C(O)COP([O-])([O-])=O", // 230,"sorbitol-6-phosphate|D-mannitol-1-phosphate|",">CH- : 4 | -OH : 5 | -CH2- : 2 | -OPO3 : 1 | Origin : 1"
		"-394.55","CC1C(O)C(O)C(O)C(O1)OP([O-])([O-])=O", // 231,"L-fuculose-1-phosphate|CHEBI:16647","-OPO3 : 1 | ring -CH< : 5 | ring -O- : 1 | -OH : 3 | -CH3 : 1 | Origin : 1"
		"-393.5","OCC(=O)C(O)C(O)COP([O-])([O-])=O", // 232,"ribulose-5-phosphate|xylulose-5-phosphate|xu5p-D|",">CH- : 2 | >C=O : 1 | -OH : 3 | -CH2- : 2 | -OPO3 : 1 | Origin : 1"
		"-479.18","C1(=C(C(=O)[O-])N(C(NC1(=O))=O)C2(C(O)C(C(COP([O-])([O-])=O)O2)O))", // 233,"orotidine-5'-phosphate|CHEBI:15842","ring -CH< : 4 | ring -O- : 1 | ring =C< : 1 | ring >C=O : 2 | -OH : 2 | ring =CH- : 1 | -COO : 1 | -CH2- : 1 | -OPO3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 2"
		"-370.27","CC(O)(CCOP([O-])([O-])=O)CC([O-])=O", // 234,"5-phosphomevalonate|CHEBI:17436","-CH2- : 3 | >C< : 1 | -COO : 1 | -CH3 : 1 | -OH : 1 | -OPO3 : 1 | Origin : 1"
		"-392.04","OC(COP([O-])([O-])=O)C(O)C(O)C=O", // 235,"D-arabinose-5-phosphate|",">CH- : 3 | -OH : 3 | -CH2- : 1 | -CH=O : 1 | -OPO3 : 1 | Origin : 1"
		"-430.78","OC1OC(COP([O-])([O-])=O)C(O)C(O)C1O", // 236,"D-galactose-6-phosphate|D-glucose-6-phosphate|D-mannose-6-phosphate|","-OPO3 : 1 | -CH2- : 1 | ring -CH< : 5 | ring -O- : 1 | -OH : 4 | Origin : 1"
		"-393.69","CC(O)C(O)C(O)C(=O)COP([O-])([O-])=O", // 238,"L-rhamnulose-1-phosphate|","-OPO3 : 1 | -CH2- : 1 | >C=O : 1 | >CH- : 3 | -OH : 3 | -CH3 : 1 | Origin : 1"
		"-85.85","CC([N+])CC([N+])CC([O-])=O", // 242,"L-erythro-3,5-diaminohexanoate|","-CH2- : 2 | >CH- : 2 | -NH3 : 2 | -CH3 : 1 | -COO : 1 | Origin : 1"
		"-378.44","OCC1OC(OC2C(O)C(O)C(O)C(O)C2O)C(O)C(O)C1O", // 245,"1-alpha-D-galactosyl-myo-inositol|","ring -CH< : 11 | -O- : 1 | -OH : 9 | ring -O- : 1 | -CH2- : 1 | Origin : 1"

// only 2 heteroaromatic rings by jankowski
// no OCCC by jankowski
//		"-529.37","Cc1cc2nc3c(nc(=O)[nH]c3=O)n(CC(O)C(O)C(O)COP([O-])(=O)OP([O-])(=O)OCC3OC(C(O)C3O)n3cnc4c(N)ncnc34)c2cc1C", // 250,"flavin-adenine-dinucleotide|fad|","ring -N< : 2 | two fused rings >C= : 4 | aromatic ring fused to nonaromatic ring >C=  : 2 | -CH2- : 3 | ring =N- : 5 | aromatic ring =CH- : 2 | >CH- : 3 | ring >C=O : 2 | aromatic ring =C< : 2 | -OH : 5 | -CH3 : 2 | -OPO3- : 1 | -OPO2- : 1 | ring -CH< : 4 | ring -O- : 1 | ring =CH- : 2 | ring =C< : 1 | -NH2 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCN : 1 | OCNC : 1 | NCCN : 1 | CCNC : 1 | heteroaromatic : 2"

		"-694.61","OCC1OC(OCC2OC(OCC3OC(OC4(CO)OC(CO)C(O)C4O)C(O)C(O)C3O)C(O)C(O)C2O)C(O)C(O)C1O", // 257,"stachyose|","ring >C< : 1 | -O- : 3 | ring -O- : 4 | ring -CH< : 18 | -CH2- : 5 | -OH : 14 | Origin : 1"

// wrong SMILES -> no "ring -NH-"
//		"-478.33","C1(=NC(=C2(C(=N1)N(C=N2)C3(OC(C(C3O)OP([O-])([O-])=O)COP([O-])(=O)[O-])))N)", // 259,"adenosine-3',5'-biphosphate|","dbl sgl ring =N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 2 | -OH : 1 | ring =C< : 1 | -OPO3 : 2 | -CH2- : 1 | ring -NH- : 1 | -NH2 : 1 | Origin : 1 | CCCN : 1 | CCNC : 3"

// only 1 heteroaromatic ring by jankowski
//		"-151.93","OCC1OC(C(O)C1O)n1cnc2c1[nH]c(=O)[nH]c2=O", // 266,"xanthosine|","two fused rings >C= : 2 | ring -N< : 1 | ring -CH< : 4 | ring =CH- : 1 | ring =N- : 1 | ring >C=O : 2 | ring -O- : 1 | -OH : 3 | -CH2- : 1 | amide : 2 | ring -NH- : 2 | Origin : 1 | OCCC : 1 | heteroaromatic : 1"

		"-337.52","C12(N(CC(O)C(O)C(O)COP([O-])(=O)[O-])c3(c(N=C2(C(=O)NC(=O)N=1))cc(C)c(C)c3))", // 270,"FMN|fmnox|","ring -N< : 1 | aromatic ring fused to nonaromatic ring >C=  : 2 | two fused rings >C= : 2 | -CH2- : 2 | aromatic ring =CH- : 2 | ring =N- : 2 | >CH- : 3 | aromatic ring =C< : 2 | ring >C=O : 2 | -OH : 3 | -CH3 : 2 | -OPO3 : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCN : 1 | OCNC : 1 | NCCN : 1 | CCNC : 1"
		"-110.12","[H]C1(CCC(=O)N1)C([O-])=O", // 280,"5-oxo-D-proline|","ring -CH< : 1 | ring -CH2- : 2 | -COO : 1 | ring >C=O : 1 | amide : 1 | ring -NH- : 1 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"-673.08","Nc1nc2n(cnc2c(=O)[nH]1)C1OC(COP([O-])(=O)OP([O-])(=O)OC2OC(CO)C(O)C(O)C2O)C(O)C1O", // 281,"GDPgalactose|GDPmannose|GDPglucose|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 9 | ring =CH- : 1 | ring =N- : 2 | ring -O- : 2 | ring >C=O : 1 | ring =C< : 1 | -OH : 6 | -NH2 : 1 | -CH2- : 2 | -OPO3- : 1 | -OPO2- : 1 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

		"-863.85","CC(C(O)=O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 288,"methylmalonyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 6 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 2 | -CH3 : 3 | >C=O : 3 | -COO : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-202.8","CC(C([O-])=O)C(C)(O)C([O-])=O", // 295,"R,S-2,3-dimethylmalate|2,3-dimethylmalate|",">C< : 1 | >CH- : 1 | -COO : 2 | -CH3 : 2 | -OH : 1 | Origin : 1"
		"-166.28","OC1CC(=CC(=O)C1O)C([O-])=O", // 297,"3-dehydroshikimate|","ring =C< : 1 | ring -CH2- : 1 | ring =CH- : 1 | -COO : 1 | ring -CH< : 2 | ring >C=O : 1 | -OH : 2 | Origin : 1 | OCCC : 2"
		"-164.73","OC1C=C(CC(=O)C1O)C([O-])=O", // 298,"5-dehydroshikimate|","ring =C< : 1 | ring -CH2- : 1 | ring =CH- : 1 | -COO : 1 | ring >C=O : 1 | ring -CH< : 2 | -OH : 2 | Origin : 1 | OCCC : 1"

// wrong SMILES -> no -COO
//		"-279.89","COc1cc(cc(OC)c1O)C=CC(=O)OC1OC(CO)C(O)C(O)C1O", // 306,"1-sinapoyl-D-glucose|","aromatic ring =C< : 4 | -O- : 3 | aromatic ring =CH- : 2 | ring -CH< : 5 | -CH3 : 2 | ring -O- : 1 | =CH- : 2 | -OH : 4 | -CH2- : 1 | -COO : 1 | Origin : 1 | OCCC : 1 | CCCC : 1"

// no heteroaromatic ring by jankowski
//		"-39.72","CC(O)C(O)C1=Nc2c(NC1)[nH]c(N)nc2=O", // 308,"7,8-dihydrobiopterin|","two fused rings >C= : 2 | ring -NH- : 2 | ring =N- : 2 | ring >C=O : 1 | ring -CH2- : 1 | ring =C< : 2 | -NH2 : 1 | >CH- : 2 | -OH : 2 | -CH3 : 1 | amide : 1 | Origin : 1 | OCCC : 1 | CCNC : 2"

		"-380.65","[N+]C(CC(=O)OP([O-])([O-])=O)C([O-])=O", // 315,"4-phospho-L-aspartate|L-4-aspartyl-phosphate|","-CO-OPO3- : 1 | -CH2- : 1 | >CH- : 1 | -COO : 1 | -NH3 : 1 | Origin : 1"

// is heteroaromatic ring by jankowski
//		"-125.24","NC(=O)C1=CN(C=CC1)C1OC(CO)C(O)C1O", // 317,"nicotinamide-riboside|","ring -CH< : 4 | dbl sgl ring =N< : 1 | ring -O- : 1 | ring =CH- : 4 | -OH : 3 | ring =C< : 1 | -CH2- : 1 | >C=O : 1 | amide : 1 | -NH2 : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 1"

		"-215.78","[H]C1(OC(=O)C(O)C1O)C(O)CO", // 324,"D-galactono-1,4-lactone|","ring -CH< : 3 | ring -O-CO- : 1 | >CH- : 1 | -OH : 4 | -CH2- : 1 | Origin : 1"
		"-291.21","CC(O)C(=O)SCC(NC(=O)CCC([N+])C(O)=O)C(=O)NCC([O-])=O", // 328,"R-S-lactoylglutathione|",">CH- : 3 | -CH2- : 4 | >C=O : 3 | -COO : 2 | -CH3 : 1 | -OH : 1 | -NH3 : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1"
		"-80.65","CC(=O)NC(Cc1ccccc1)C([O-])=O", // 332,"N-acetyl-L-phenylalanine|",">CH- : 1 | -CH2- : 1 | -COO : 1 | aromatic ring =C< : 1 | >C=O : 1 | aromatic ring =CH- : 5 | -CH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"

// wrong SMILES -> no ring
//		"-248.86","[H]C(O)(CO)C([H])(O)C([H])(O)C([H])(NC(=O)CO)C=O", // 333,"N-glycolyl-D-mannosamine|","ring -CH< : 5 | -OH : 5 | ring -O- : 1 | >C=O : 1 | -CH2- : 2 | amide : 1 | -NH- : 1 | Origin : 1"

// wrong SMILES -> no -NH3 possible
//		"-165.94","CNC(CC(O)=O)C([O-])=O", // 336,"L-threo-3-methylaspartate|",">CH- : 2 | -COO : 2 | -CH3 : 1 | -NH3 : 1 | Origin : 1"

		"-153.08","OC(=O)CC1NC(=O)NC1=O", // 340,"L-5-carboxymethylhydantoin|","ring -CH< : 1 | ring >C=O : 2 | -CH2- : 1 | -COO : 1 | amide : 2 | ring -NH- : 2 | Origin : 1"
		"-36.26","[H]C12CCC3([H])C4([H])CCC(=O)C4(C)CCC3([H])C1(C)CCC(=O)C2", // 342,"5alpha-androstane-3,17-dione|5beta-androstane-3,17-dione|","two fused rings >CH- : 4 | ring -CH2- : 9 | two fused rings >C< : 2 | -CH3 : 2 | ring >C=O : 2 | Origin : 1"
		"-408.89","OC1C(COP([O-])([O-])=O)OC(C1O)n1cnc2c(NC(CC(O)=O)C(O)=O)ncnc12", // 343,"adenylosuccinate|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 4 | ring =CH- : 2 | ring =N- : 3 | ring -O- : 1 | ring =C< : 1 | -OH : 2 | -NH- : 1 | -CH2- : 2 | >CH- : 1 | -OPO3 : 1 | -COO : 2 | Origin : 1 | heteroaromatic : 2"

// wrong SMILES -> no -OPO3 by jankowski
//		"0.09","NC(=O)c1ncn(C2OC(COP([O-])([O-])=O)C(O)C2O)c1N", // 349,"5-amino-4-imidazolecarboxamide|","ring =C< : 2 | >C=O : 1 | ring =N- : 1 | ring -NH- : 1 | -NH2 : 2 | ring =CH- : 1 | amide : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 1"

		"-183.72","OCC(O)CC(=O)C([O-])=O", // 351,"2-dehydro-3-deoxy-L-pentonate|","-COO : 1 | >C=O : 1 | -CH2- : 2 | >CH- : 1 | -OH : 2 | Origin : 1 | OCCO : 1"
		"-425.13","CC(=O)NC1C(O)C(O)C(CO)OC1OP([O-])([O-])=O", // 355,"N-acetyl-alpha-D-glucosamine-1-phosphate|N-acetyl-D-glucosamine-1-phosphate|","ring -CH< : 5 | ring -O- : 1 | -OPO3 : 1 | -OH : 3 | >C=O : 1 | -CH3 : 1 | -CH2- : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-425.13","CC(=O)NC1(C(O)OC(COP([O-])([O-])=O)C(O)C1(O))", // 356,"N-acetyl-D-mannosamine-6-phosphate|N-acetyl-D-glucosamine-6-phosphate|","ring -CH< : 5 | -OH : 3 | ring -O- : 1 | >C=O : 1 | -CH3 : 1 | -CH2- : 1 | -OPO3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"

// no heteroaromatic ring by jankowski
//		"-651.64","CC1OC(OP([O-])(=O)OP([O-])(=O)OCC2OC(CC2O)n2cc(C)c(=O)[nH]c2=O)C(O)C(O)C1N", // 357,"dTDP-4-amino-4,6-dideoxy-D-glucose|","ring -CH< : 8 | ring -O- : 2 | ring -CH2- : 1 | ring >C=O : 2 | ring =CH- : 1 | ring =C< : 1 | -CH2- : 1 | -OH : 3 | -CH3 : 2 | -OPO3- : 1 | -OPO2- : 1 | -NH3 : 1 | amide : 2 | ring -N< : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1"

		"-49.38","[H]C12CCC3([H])C4([H])CCC(O)C4(C)CCC3([H])C1(C)CCC(O)C2", // 363,"5alpha-androstane-3alpha,17beta-diol|5alpha-androstane-3beta,17alpha-diol|","two fused rings >CH- : 4 | ring -CH2- : 9 | two fused rings >C< : 2 | -CH3 : 2 | ring -CH< : 2 | -OH : 2 | Origin : 1"

// wrong SMILES -> same as 363
//		"-42.82","[H]C12CCC3([H])C4([H])CCC(=O)C4(C)CCC3([H])C1(C)CCC(O)C2", // 364,"5alpha-androstane-3alpha-ol-17-one|5beta-androstane-3alpha-ol-17-one|","two fused rings >CH- : 4 | ring -CH2- : 9 | two fused rings >C< : 2 | -CH3 : 2 | ring >C=O : 1 | ring -CH< : 1 | -OH : 1 | Origin : 1"

		"-294.26","[O-]C(=O)CCC(=O)NC(CCCC(=O)C([O-])=O)C([O-])=O", // 365,"N-succinyl-2-L-amino-6-oxoheptanedioate|",">CH- : 1 | -CH2- : 5 | -COO : 3 | >C=O : 2 | amide : 1 | -NH- : 1 | Origin : 1 | OCCO : 1"
		"-432.64","OC(COP([O-])([O-])=O)C(O)CC(=O)C([O-])=O", // 366,"6-phospho-2-dehydro-3-deoxy-D-gluconate|2-dehydro-3-deoxy-D-galactonate-6-phosphate|","-OPO3 : 1 | -CH2- : 2 | >CH- : 2 | -OH : 2 | >C=O : 1 | -COO : 1 | Origin : 1 | OCCO : 1"
		"-269.49","[N+]C(CCCC(NC(=O)CCC(O)=O)C(O)=O)C([O-])=O", // 367,"N-succinyl-L-2,6-diaminoheptanedioate|",">CH- : 2 | -CH2- : 5 | -COO : 3 | >C=O : 1 | -NH3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-817.39","CCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 372,"3S-3-hydroxyhexanoyl-CoA|S-3-hydroxyhexanoyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 3 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 9 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 2 | -CH3 : 3 | >C=O : 3 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-809.37","CCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 373,"3-oxohexanoyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 9 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 4 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-758.54","CCCC=CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 374,"cis-hex-2-enoyl-CoA|trans-hex-2-enoyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 8 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 3 | =CH- : 2 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 2"
		"-142.05","[H]C12CCC3([H])C4([H])CCC(O)(C(=O)CO)C4(C)CC(=O)C3([H])C1(C)CCC(=O)C2", // 376,"5alpha-pregnane-17alpha,21-diol-3,11,20-trione|5beta-pregnane-17alpha,21-diol-3,11,20-trione|","two fused rings >CH- : 4 | two fused rings >C< : 2 | ring -CH2- : 8 | ring >C< : 1 | -CH3 : 2 | ring >C=O : 2 | >C=O : 1 | -OH : 2 | -CH2- : 1 | Origin : 1"

// only 1 heteroaromatic ring by jankowski
//		"-61.24","OCC1OC(CC1O)n1cnc2c1nc[nH]c2=O", // 380,"2'-deoxyinosine|","ring -N< : 1 | two fused rings >C= : 2 | ring -CH< : 3 | ring =CH- : 2 | ring =N- : 2 | ring -O- : 1 | ring -CH2- : 1 | ring >C=O : 1 | -CH2- : 1 | -OH : 2 | amide : 1 | ring -NH- : 1 | Origin : 1 | OCCC : 1 | CCNC : 1 | heteroaromatic : 1"

// only 1 amide by jankowski : WHO IS RIGHT ?
// --> are lactams not counted as amides ?
//		"-120.89","[H]C12SC(C)(C)C(N1C(=O)C2NC(=O)C(N)c1ccc(O)cc1)C([O-])=O", // 382,"amoxicillin|","two fused rings -N< : 1 | two fused rings >CH- : 1 | ring -CH< : 2 | ring >C=O : 1 | ring -S- : 1 | ring >C< : 1 | -COO : 1 | -CH3 : 2 | >C=O : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | >CH- : 1 | -NH2 : 1 | -OH : 1 | amide : 1 | -NH- : 1 | Origin : 1"

		"-345.38","C(C2(C(C(C(N1(C=CC(N)=NC1=O))O2)O)OP([O-])([O-])=O))O", // 383,"cytidine-3'-monophosphate|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -OPO3 : 1 | -CH2- : 1 | ring =C< : 1 | -NH2 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-217.61","C1(CO)(O)OCC(O)C(O)C1(O)", // 392,"D-psicose|L-sorbose|D-tagatose|","ring >C< : 1 | ring -CH< : 3 | ring -O- : 1 | -CH2- : 1 | -OH : 5 | ring -CH2- : 1 | Origin : 1"
		"-152.6","OCC(O)C(O)CO", // 402,"L-threitol|","-OH : 4 | -CH2- : 2 | >CH- : 2 | Origin : 1"

// ERROR : problem about maximal ring size : should be upper bounded otherwise "-O-" is never matched
//		"-950.88","OCC1OC2OC3C(CO)OC(OC4C(CO)OC(OC5C(CO)OC(OC6C(CO)OC(OC7C(CO)OC(OC1C(O)C2O)C(O)C7O)C(O)C6O)C(O)C5O)C(O)C4O)C(O)C3O", // 415,"cyclomaltohexaose|","ring -CH< : 30 | -O- : 6 | ring -O- : 6 | -CH2- : 6 | -OH : 18 | Origin : 1"
//		"-1267.84","OCC1OC2OC3C(CO)OC(OC4C(CO)OC(OC5C(CO)OC(OC6C(CO)OC(OC7C(CO)OC(OC8C(CO)OC(OC9C(CO)OC(OC1C(O)C2O)C(O)C9O)C(O)C8O)C(O)C7O)C(O)C6O)C(O)C5O)C(O)C4O)C(O)C3O", // 420,"cyclomaltooctaose|","ring -CH< : 40 | -O- : 8 | ring -O- : 8 | -CH2- : 8 | -OH : 24 | Origin : 1"
//		"-532.92","OCC1OC(OC2C(O)C(CO)OC(OC3C(O)C(O)OC(CO)C3O)C2O)C(O)C(O)C1O", // 423,"laminaritriose|beta-laminaritriose|","-OH : 11 | ring -CH< : 10 | -O- : 2 | ring -O- : 2 | -CH2- : 3 | -CH=O : 1 | >CH- : 4 | Origin : 1"

		"-336.65","[N+]C(COP([O-])([O-])=O)C([O-])=O", // 436,"L-O-phosphoserine|",">CH- : 1 | -COO : 1 | -NH3 : 1 | -CH2- : 1 | -OPO3 : 1 | Origin : 1"
		"-673.85","Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C1O", // 439,"ATP|atp|","ring =CH- : 2 | ring =N- : 3 | two fused rings >C= : 2 | ring =C< : 1 | -NH2 : 1 | ring -N< : 1 | ring -CH< : 4 | -OH : 2 | ring -O- : 1 | -OPO2- : 2 | -OPO3 : 1 | -CH2- : 1 | Origin : 1 | heteroaromatic : 2"
		"-465.85","Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C1O", // 440,"ADP|adp|","ring =CH- : 2 | two fused rings >C= : 2 | ring =C< : 1 | -NH2 : 1 | ring =N- : 3 | ring -N< : 1 | ring -CH< : 4 | -OH : 2 | ring -O- : 1 | -CH2- : 1 | -OPO2- : 1 | -OPO3 : 1 | Origin : 1 | heteroaromatic : 2"
		"-218.01","Nc1ncnc2n(cnc12)C1CC(O)C(COP([O-])([O-])=O)O1", // 453,"2'-deoxyadenosine-5'-monophosphate|dAMP|","ring =C< : 1 | ring =CH- : 2 | -NH2 : 1 | ring =N- : 3 | two fused rings >C= : 2 | ring -CH< : 3 | ring -CH2- : 1 | ring -O- : 1 | -CH2- : 1 | -OH : 1 | -OPO3 : 1 | ring -N< : 1 | Origin : 1 | heteroaromatic : 2"
		"-367.85","Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OS([O-])(=O)=O)C(O)C1O", // 475,"adenylyl-sulfate|aps|","ring =CH- : 2 | two fused rings >C= : 2 | ring =N- : 3 | ring =C< : 1 | -NH2 : 1 | ring -N< : 1 | ring -CH< : 4 | ring -O- : 1 | -OH : 2 | -CH2- : 1 | -OPO2- : 1 | -SO4 : 1 | Origin : 1 | heteroaromatic : 2"
		"-214.49","OC1CC(O)(CC(=O)C1O)C([O-])=O", // 479,"3-dehydroquinate|5-dehydroquinate|","ring -CH2- : 2 | ring -CH< : 2 | -OH : 3 | ring >C=O : 1 | ring >C< : 1 | -COO : 1 | Origin : 1"

// no heteroaromatic ring by jankowski
//		"-825.37","CCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 481,"3R-3-hydroxyhexanoyl-CoA|","dbl sgl ring =N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 2 | -OH : 3 | ring =C< : 1 | -OPO3 : 1 | -CH2- : 9 | ring -NH- : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 2 | -CH3 : 3 | >C=O : 3 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | CCCN : 1 | CCNC : 3"

		"-179.54","[H]C(O)(CO)C([H])(O)C([H])(O)C=O", // 488,"D-lyxose|","-OH : 4 | -CH=O : 1 | >CH- : 3 | -CH2- : 1 | Origin : 1"
		"-589.26","OCC1OC(OC2OC(COP([O-])([O-])=O)C(O)C(O)C2O)C(O)C(O)C1O", // 493,"alpha,alpha-trehalose-6-phosphate|","ring -CH< : 10 | -CH2- : 2 | -OPO3 : 1 | -OH : 7 | ring -O- : 2 | -O- : 1 | Origin : 1"
		"-820.63","CC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 498,"3-hydroxybutanoyl-CoA|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 3 | ring =C< : 1 | -CH2- : 7 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 2 | -CH3 : 3 | >C=O : 3 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"
		"-782.59","CCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12", // 501,"butyryl-CoA|btcoa|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 8 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 3 | >C=O : 3 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | heteroaromatic : 2"

// no heteroaromatic ring by jankowski
//		"-113.97","Nc1nc2NCC(CNc3ccc(cc3)C(=O)NC(CCC([O-])=O)C([O-])=O)Nc2c(=O)[nH]1", // 511,"5,6,7,8-tetrahydrofolate|thf|","-CH2- : 3 | -NH- : 2 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | >C=O : 1 | >CH- : 1 | -COO : 2 | ring -CH2- : 1 | two fused rings >C= : 2 | ring >C=O : 1 | ring =C< : 1 | ring -NH- : 3 | ring =N- : 1 | ring -CH< : 1 | -NH2 : 1 | amide : 2 | Origin : 1 | OCCC : 2 | CCNC : 1"
//		"-93.65","Nc1nc2NCC3CN(CN3c2c(=O)[nH]1)c1ccc(cc1)C(=O)NC(CCC([O-])=O)C([O-])=O", // 513,"5,10-methylenetetrahydrofolate|","two fused rings -N< : 1 | two fused rings >C= : 2 | two fused rings >CH- : 1 | ring -CH2- : 3 | ring >C=O : 1 | ring -N< : 1 | ring -NH- : 2 | ring =N- : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | >C=O : 1 | >CH- : 1 | -CH2- : 2 | -COO : 2 | ring =C< : 1 | -NH2 : 1 | amide : 2 | -NH- : 1 | Origin : 1 | OCCC : 2 | CCNC : 1"
//		"-144.77","[H]C(=O)N(CC1CNc2nc(N)[nH]c(=O)c2N1)c1ccc(cc1)C(=O)NC(CCC([O-])=O)C([O-])=O", // 514,"10-formyltetrahydrofolate|10fthf|","two fused rings >C= : 2 | ring -NH- : 3 | ring >C=O : 1 | ring =N- : 1 | ring -CH< : 1 | ring -CH2- : 1 | ring =C< : 1 | -CH2- : 3 | -NH2 : 1 | aromatic ring =C< : 2 | -CH=O : 1 | aromatic ring =CH- : 4 | >C=O : 1 | >CH- : 1 | -COO : 2 | amide : 3 | -N< : 1 | -NH- : 1 | Origin : 1 | OCCC : 2 | CCNC : 1"
//		"-106.95","[H]C(=N)N1C(CNc2ccc(cc2)C(=O)NC(CCC(O)=O)C(O)=O)CNc2nc(N)[nH]c(=O)c12", // 516,"5-formiminotetrahydrofolate|","two fused rings >C= : 2 | ring -N< : 1 | ring >C=O : 1 | ring -CH< : 1 | =CH- : 1 | ring -NH- : 2 | ring =N- : 1 | ring -CH2- : 1 | -CH2- : 3 | =NH : 1 | ring =C< : 1 | -NH- : 2 | -NH2 : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | >C=O : 1 | >CH- : 1 | -COO : 2 | amide : 2 | Origin : 1 | OCCC : 2 | CCNC : 1"
//		"-109.12","Nc1nc2NCC(CNc3ccc(cc3)C(=O)NC(CCC(O)=O)C(O)=O)=Nc2c(=O)[nH]1", // 517,"7,8-dihydrofolate|","two fused rings >C= : 2 | ring -NH- : 2 | ring =N- : 2 | ring >C=O : 1 | ring -CH2- : 1 | ring =C< : 2 | -NH2 : 1 | -CH2- : 3 | -NH- : 2 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | >C=O : 1 | >CH- : 1 | -COO : 2 | amide : 2 | Origin : 1 | OCCC : 2 | CCNC : 2"

// no amide by jankowski : WHO IS RIGHT ?
// --> are lactams not counted as amides ?
//		"-85.89","[H]C1(N)C(=O)N2C1([H])SC(C)(C)C2([H])C([O-])=O", // 521,"6-aminopenicillanate|","two fused rings -N< : 1 | two fused rings >CH- : 1 | ring -CH< : 2 | ring >C=O : 1 | ring -S- : 1 | ring >C< : 1 | -COO : 1 | -NH2 : 1 | -CH3 : 2 | Origin : 1"

		"-376.76","OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O", // 535,"cellobiose|beta-lactose|beta-maltose|","ring -CH< : 10 | -O- : 1 | ring -O- : 2 | -CH2- : 2 | -OH : 8 | Origin : 1"
		"-430.78","OCC1OC(OP([O-])([O-])=O)C(O)C(O)C1O", // 547,"glucose-1-phosphate|galactose-1-phosphate|mannose-1-phosphate|","ring -CH< : 5 | ring -O- : 1 | -OPO3 : 1 | -OH : 4 | -CH2- : 1 | Origin : 1"
		"-764.06","CC(C)(COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCSC(=O)C=C", // 560,"acrylyl-coa|","ring -N< : 1 | ring -CH< : 4 | two fused rings >C= : 2 | ring =CH- : 2 | ring -O- : 1 | ring =N- : 3 | -OH : 2 | ring =C< : 1 | -CH2- : 6 | -OPO3 : 1 | -NH2 : 1 | -OPO3- : 1 | -OPO2- : 1 | >C< : 1 | >CH- : 1 | -CH3 : 2 | >C=O : 3 | =CH- : 1 | =CH2 : 1 | amide : 2 | -NH- : 2 | -S- : 1 | thioester : 1 | Origin : 1 | OCCC : 1 | heteroaromatic : 2"
		"-345.38","C(C2(C(C(C(N1(C(N=C(C=C1)N)=O))O2)O)O))OP([O-])([O-])=O", // 586,"cmp|","ring -CH< : 4 | ring -O- : 1 | ring >C=O : 1 | ring =CH- : 2 | -OH : 2 | ring =N- : 1 | -CH2- : 1 | ring =C< : 1 | -OPO3 : 1 | -NH2 : 1 | amide : 1 | ring -N< : 1 | Origin : 1 | OCNC : 1 | CCCN : 1"
		"-16.93","CC(Cl)CCl", // 609,"1,2-dichloropropane|","-CH3 : 1 | -CH2- : 1 | primary -Cl : 1 | >CH- : 1 | secondary -Cl : 1 | Origin : 1 | BinaryVicinalCl : 1 | DistinctVicinalCl : 1"

// wrong SMILES : no chlorines ...
//		"-25.58","[H]C1(C(C(=O)N2CCCC2OC)C2(O)C(O)C1(Oc1cc(OC)cc(OC)c21)c1ccc(OC)cc1)c1ccccc1", // 674,"2,3,4,5-tetrachlorophenol|","aromatic ring =CH- : 1 | aromatic ring =C< : 5 | -OH : 1 | tertiary -Cl : 4 | Origin : 1 | BinaryVicinalCl : 3 | DistinctVicinalCl : 3"

		"-466.4","[N+]Cc1cc(COc2ccc(CCNC(=O)CCC(NC(=O)CCC(NC(=O)CCC(C(CCC(O)=O)C(O)=O)C(O)=O)C(O)=O)C(O)=O)cc2)co1", // 683,"MFR|methanofuran|","ring =C< : 2 | ring =CH- : 2 | -CH2- : 12 | -O- : 1 | ring -O- : 1 | aromatic ring =C< : 2 | -NH3 : 1 | aromatic ring =CH- : 4 | >C=O : 3 | >CH- : 4 | -COO : 5 | amide : 3 | -NH- : 3 | Origin : 1 | heteroaromatic : 1"

// wrong SMILES : no -COO, no -OPO3-
//		"-544.24","C1(NC2(N=C(N)NC(=O)C=2(NC1C(O)C(O)CO)))", // 685,"H4MPT|","two fused rings >C= : 2 | ring -NH- : 3 | ring >C=O : 1 | ring =N- : 1 | ring -CH< : 6 | ring =C< : 1 | >CH- : 5 | -NH2 : 1 | -NH- : 1 | -CH3 : 2 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -CH2- : 5 | -OH : 5 | -O- : 1 | ring -O- : 1 | -OPO3- : 1 | -COO : 2 | amide : 1 | Origin : 1 | OCCC : 1 | CCNC : 1"

// no heteroaromatic ring by jankowski
//		"-529.04","[H]C1(Nc2c(NC1C)nc(N)[nH]c2=O)C(C)Nc1ccc(CC(O)C(O)C(O)COC2OC(COP([O-])(=O)OC(CCC(O)=O)C(O)=O)C(O)C2O)cc1", // 687,"methenyl-H4MPT|","dbl sgl ring =N< : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -CH2- : 5 | >CH- : 4 | -OH : 5 | -O- : 1 | ring -CH< : 6 | ring -O- : 1 | -OPO3- : 1 | -COO : 2 | two fused rings >C= : 2 | ring -NH- : 2 | ring =N- : 1 | two fused rings -N< : 1 | ring >C=O : 1 | ring =C< : 1 | two fused rings >CH- : 1 | ring =CH- : 1 | -CH3 : 2 | -NH2 : 1 | amide : 1 | Origin : 1 | OCCC : 1 | CCNC : 2"
//		"-653.72","[H]C(CCC([O-])=O)(NC(=O)CCC([H])(NC(=O)C(C)OP([O-])(=O)OCC(O)C(O)C(O)CN1c2cc(O)ccc2Cc2c1[nH]c(=O)[nH]c2=O)C([O-])=O)C([O-])=O", // 689,"F420|","aromatic ring =CH- : 3 | aromatic ring =C< : 1 | aromatic ring fused to nonaromatic ring >C=  : 2 | ring -N< : 1 | ring =CH- : 1 | two fused rings >C= : 2 | ring =N- : 1 | ring >C=O : 2 | -OH : 4 | -CH2- : 6 | >CH- : 6 | -OPO3- : 1 | >C=O : 2 | -CH3 : 1 | -COO : 3 | amide : 3 | ring -NH- : 1 | -NH- : 2 | Origin : 1 | OCCC : 1 | OCNC : 1 | CCCN : 1 | CCCC : 1"
//		"-531.97","[H]C1(C(C)Nc2ccc(CC(O)C(O)C(O)COC3OC(COP([O-])(=O)OC(CCC(O)=O)C(O)=O)C(O)C3O)cc2)C(C)Nc2nc(N)[nH]c(=O)c2N1C", // 691,"CH3-H4MPT|","two fused rings >C= : 2 | ring -NH- : 2 | ring =N- : 1 | ring -N< : 1 | ring >C=O : 1 | ring -CH< : 6 | ring =C< : 1 | -CH3 : 3 | -NH2 : 1 | >CH- : 5 | -NH- : 1 | aromatic ring =C< : 2 | aromatic ring =CH- : 4 | -CH2- : 5 | -OH : 5 | -O- : 1 | ring -O- : 1 | -OPO3- : 1 | -COO : 2 | amide : 1 | Origin : 1 | OCCC : 1 | CCNC : 1"

		"-353.81","CC(OP([O-])([O-])=O)C(NC(=O)CCCCCCS)C([O-])=O", // 695,"CoB|H-S-CoB|coenzyme-B|",">C=O : 1 | -CH3 : 1 | >CH- : 2 | -OPO3 : 1 | -COO : 1 | -CH2- : 6 | -SH : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-467.89","CC(OP([O-])([O-])=O)C(NC(=O)CCCCCCSSCCS([O-])(=O)=O)C([O-])=O", // 696,"CoM-S-S-CoB|",">C=O : 1 | -CH3 : 1 | >CH- : 2 | -OPO3 : 1 | -COO : 1 | -CH2- : 8 | -S-S- : 1 | -SO3 : 1 | amide : 1 | -NH- : 1 | Origin : 1"
		"-376.76","OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O", // 697,"O-alpha-D-glucopyranosyl-1->4-n|maltose|disaccharide|lactose|lcts|","ring -CH< : 10 | -O- : 1 | ring -O- : 2 | -CH2- : 2 | -OH : 8 | Origin : 1"
		"-376.76","OCC1OC(OCC2OC(O)C(O)C(O)C2O)C(O)C(O)C1O", // 698,"allolactose|","ring -CH< : 10 | -CH2- : 2 | -OH : 8 | ring -O- : 2 | -O- : 1 | Origin : 1"
		"-535.24","OCC1OC(OC2C(CO)OC(OC3C(CO)OC(O)C(O)C3O)C(O)C2O)C(O)C(O)C1O", // 708,"cellotriose|beta-cellotriose|","ring -CH< : 15 | -O- : 2 | ring -O- : 3 | -CH2- : 3 | -OH : 11 | Origin : 1"
		"-377.65","OCC1OC(OC2C(CO)OC(O)(CO)C2O)C(O)C(O)C1O", // 709,"lactulose|beta-lactulose|","ring -O- : 2 | ring -CH< : 8 | -O- : 1 | -OH : 8 | ring >C< : 1 | -CH2- : 3 | Origin : 1"

		"-570","P1(O)(=O)OP(O)(=O)OP(O)(=O)O1", // "cyclic trimetaphosphate"
//		"NOT HANDLED BY JANKOWSKI","OP(O)(=O)OP(O)(=O)OP(O)(=O)O", // "triphosphoric acid"

	"" // marks the end of the array
};

class DR_print : public MoleculeDecomposition::DecompositionReporter {
public:
	  /*!
	   * Reports a component that was matched and its according ID
	   * in the final match graph
	   * @param component the component that was matched
	   * @param matchCount the number of matches found
	   * @param machtID the ID of the component in the final match graph
	   */
	void
	reportComponent(	const MoleculeComponent & component,
						const size_t matchCount,
						const size_t matchID )
	{
		std::cout <<" + " << matchCount <<" [" <<matchID <<"] = " <<component.description <<std::endl;
	}

	  /*!
	   * Reports an interaction pattern that was matched and its
	   * matched indices
	   * @param interactionDescription the interaction that was matched
	   * @param matchCount the number of matches found
	   */
	virtual
	void
	reportInteraction(	const std::string & interactionDescription,
						const size_t matchCount )
	{
		std::cout <<" * " <<matchCount <<" " <<interactionDescription <<std::endl;
	}

	  /*!
	   * Reports the matching of a full small molecule.
	   * @param smallMolecule the small molecule matched
	   */
	virtual
	void
	reportSmallMolecule( const MoleculeComponent & smallMolecule )
	{
		std::cout <<" = " <<smallMolecule.description <<std::endl;
	}

	  /*!
	   * Reports the final match graph where each matched component ID
	   * is given as class information of the matched nodes.
	   *
	   * @param graph the matched graph
	   */
	void
	reportMatchGraph(	const sgm::Graph_Interface & graph )
	{
		std::cout <<"\n = final match graph :\n";
		print( graph );
		std::cout <<std::endl;
	}

	void
	reportMatchComplete( const bool matchComplete )
	{
		std::cout <<"\n = match complete = " <<(matchComplete?"true":"false");
		std::cout <<std::endl;
	}

	void
	reportPolyPhosphate( const size_t numOfMiddlePhosphates
									, const size_t matchID )
	{
		std::cout <<" + " << numOfMiddlePhosphates
					<<" [" <<matchID <<"] = middle phosphate '-O-PO2^{-1}-'"
					<<std::endl;
	}
};



int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"      ggl::chem::MoleculeDecomposition        \n"
				<<"==============================================\n" 
				<<std::endl;

	sgm::GM_vf2 gmVF2;
	sgm::SGM_vf2 sgmVF2;

	for (size_t i = 0; !(SMILES[i].empty()); i++) {

		  // even entries hold energy value
		if (i%2 == 0) {
			continue;
		}

		std::cout	<<"\n####### next SMILES ########\n\n"
					<<" SMILES = " <<SMILES[i]
					<<std::endl;

		std::cout <<"\n-->SMILESparser::parseSMILES( SMILES )" <<std::endl;

		std::pair<Molecule, int> ret = SMILESparser::parseSMILES(SMILES[i]);

		if (ret.second < 0 )
		{
			std::cout <<"\n Molecule GRAPH to be evaluated : " <<std::endl;
			std::cout <<ret.first;
			std::cout <<"\n" <<std::endl;
		} else {
			std::cout <<"\n PARSING ERROR : at input position " <<ret.second
					<<" = '" <<SMILES[i].at(ret.second) <<"' "
					<<" within \n'"
					<<SMILES[i].substr(std::max(0,((int)ret.second)-10),20)
					<<"'\n"
					<<std::endl;
		}

		std::cout << " MOLECULE parsed = \n"
					<<Molecule_Graph(ret.first)
					<<std::endl;

		  // add missing protons
		Molecule molAllProtons;
		MoleculeUtil::copy( ret.first, molAllProtons );
		MoleculeUtil::fillProtons( molAllProtons );

		std::cout << " MOLECULE corrected = \n"
					<<Molecule_Graph(molAllProtons)
					<<std::endl;

		DR_print reporter;
		MoleculeDecomposition energyEval(gmVF2, sgmVF2, reporter);

//std::cerr <<"\n : "<<SMILES[i] <<"\n";
		std::cout <<"\n-->MoleculeDecomposition::getEnergy( MOLECULE corrected ) :\n" <<std::endl;
		double energy = energyEval.getEnergy( molAllProtons );
		std::cout <<"\n = final energy = "<<energy <<" kcal/mol"<<std::endl;
		std::cout <<"\n = expected was = "<<SMILES[i-1] <<" kcal/mol"<<std::endl;
	}



	std::cout	<<"\n"
				<<"===============  END TEST  ===================\n"
				<<std::endl;
	return 0;
}
