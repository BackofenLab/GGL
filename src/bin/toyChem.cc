

#include "toyChemUtil.hh"

#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cstring>



#include "biu/OptionParser.hh"


#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"
#include "sgm/PA_OrderCheck.hh"


#include "ggl/MR_ApplyRule.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/GS_stream.hh"

#include <ggl/chem/SMILESparser.hh>
#include "ggl/chem/GS_SMILES.hh"
#include "ggl/chem/SMILESwriter.hh"

#include "ggl/chem/Reaction.hh"
#include "ggl/chem/MR_Reactions.hh"
#include "ggl/chem/RRC_ArrheniusLaw.hh"
#include "ggl/chem/EC_MoleculeDecomposition.hh"
#include "RRC_QuantumMechanics.hh"

#include "ggl/chem/AP_NSPDK.hh"
#include "ggl/chem/AP_disabled.hh"

#include "version.hh"

//////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
#if HAVE_UNORDERED_MAP > 0
#include <unordered_map>
#define USED_HASH_MAP_DESCRIPTION "std::unordered_map"

#elif HAVE_TR1_UNORDERED_MAP > 0
#include <tr1/unordered_map>
#define USED_HASH_MAP_DESCRIPTION "std::tr1::unordered_map"

#elif HAVE_GNU_HASH_MAP > 0
#include <ext/hash_map>
#define USED_HASH_MAP_DESCRIPTION "__gnu_cxx::hash_map"

#else
#include <map>
#define USED_HASH_MAP_DESCRIPTION "std::map"

#endif


//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments( biu::OptionMap & allowedArgs, std::string &infoText );  


//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraph( std::ostream& out, ggl::chem::Molecule & m );

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphs( std::ostream& out, SMILES_container & smiles );

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphsGML( std::ostream& out, SMILES_container & smiles, const bool withSpaces );

//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample();

//////////////////////////////////////////////////////////////////////////

void
giveGraphGMLExample();

//////////////////////////////////////////////////////////////////////////

CompareStringPointer compareStringPointer;

//////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv ) {
	using namespace std;
	using namespace ggl;
	using namespace chem;
	
	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////
	int exitValue = 0;
	
	enum InfoMode {OUT_SILENT, OUT_NORMAL, OUT_VERBOSE};
	InfoMode infoMode = OUT_NORMAL;
	
	enum PrintMode {PRINT_SMILES, PRINT_GRAPHS, PRINT_GML, PRINT_REACTIONS, PRINT_REACTION_NETWORK};
	PrintMode printMode = PRINT_SMILES;
	
	std::ostream* out = &std::cout;
	std::ofstream* outFile = NULL;
	
	ReactionRateCalculation*  rateCalc = NULL;
	EC_MoleculeDecomposition* ec_molDec = NULL;
	
	AromaticityPerception * aromaticityPrediction = NULL;

	size_t iterations = 0;
	 // rule pattern for each number of connected components
	RulePatternMap rulePattern;
	  // the supported molecule groups that can be abbreviated in rule and
	  // molecule input
	ggl::chem::GroupMap moleculeGroups;
	
	
	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////
	
	biu::OptionMap allowedArgs;	//< map of allowed arguments
	std::string infoText;		//< info string of the program
	initAllowedArguments(allowedArgs,infoText);	// init
	
		// parse programm arguments	
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
													argv, infoText);
		// arguments parseable and all mandatory arguments given
	if (opts.argExist("help")) {
		opts.coutUsage();
		return 0;
	}
	if (opts.argExist("ruleExample")) {
		giveRuleGMLExample();
		return 0;
	}
	if (opts.argExist("graphExample")) {
		giveGraphGMLExample();
		return 0;
	}
	if (opts.argExist("version")) {
		giveVersion();
		std::cout	<<"\n used SMILES_container : " <<USED_HASH_MAP_DESCRIPTION <<"\n"
					<<std::endl;
		return 0;
	}
	if (opts.noErrors()) {
		  // rules are officially not mandatory.. we have to request manually
		if (!opts.argExist("rules") || opts.getStrVal("rules").size() == 0 ) {
			std::cerr <<"\n\n\tERROR : needed argument not given : 'rules' check usage\n\n";
			return -1;
		}
		  // rules are officially not mandatory.. we have to request manually
		if (	(!opts.argExist("smiles") || opts.getStrVal("smiles").size()==0) 
			&&	(!opts.argExist("mols") || opts.getStrVal("mols").size()==0) ) 
		{
			std::cerr <<"\n\n\tERROR : needed argument not given :"
				<<" neither 'smiles' nor 'mols' is present, check usage\n\n";
			return -1;
		}
		if (opts.getBoolVal("v")) {
			infoMode = OUT_VERBOSE;
		}
	} else {
		return -1;
	}

	SMILES_container c1;
	SMILES_container c2;
	SMILES_container& targetSmiles = c1;
	SMILES_container& producedSmiles = c2;
	
	try { 


		//////////////////////////////////////////////////////////////
		// set up iteration parameters
		//////////////////////////////////////////////////////////////
		
		if (opts.getIntVal("iter") < 0) {
			std::ostringstream oss;
			oss	<<"number of rule application iterations (" 
			<<opts.getIntVal("iter") <<") has to be at least 0";
			throw ArgException(oss.str());
		}
		iterations = (size_t)opts.getIntVal("iter");
		
		const bool allowAllIntra = opts.getBoolVal("allowAllIntra");
		
		//////////////////////////////////////////////////////////////
		// check for multiple use of STDIN as input parameter value
		//////////////////////////////////////////////////////////////
		
		{
			int conflicts = 1;
			conflicts *= (opts.getStrVal("rules").compare("STDIN")==0?2:1);
			conflicts *= (opts.getStrVal("smiles").compare("STDIN")==0?3:1);
			conflicts *= (opts.getStrVal("rules").compare("STDIN")==0?5:1);
			
			if (conflicts > 5) {
				std::ostringstream oss;
				oss <<"cannot read ";
				switch (conflicts) {
					case 6: oss <<"RULES and SMILES"; break;
					case 10: oss <<"RULES and MOLECULES"; break;
					case 15: oss <<"SMILES and MOLECULES"; break;
					case 30: oss <<"RULES, SMILES and MOLECULES"; break;
				}
				oss <<" from STDIN!";
				throw ArgException(oss.str());
			}
		}
		
		//////////////////////////////////////////////////////////////
		// set up output streams
		//////////////////////////////////////////////////////////////
		
		  // set output stream
		if (opts.getStrVal("out").size() == 0) {
			throw ArgException("no output file given");
		} else if (opts.getStrVal("out").compare("STDOUT") != 0) {
			outFile = new std::ofstream(	opts.getStrVal("out").c_str()
											, std::ofstream::out );
			if (!outFile->is_open()) {
				std::ostringstream oss;
				oss	<<"cannot open output file '" <<opts.getStrVal("out") <<"'";
				throw ArgException(oss.str());
			}
			out = outFile;
		}

		//////////////////////////////////////////////////////////////
		// setup output details
		//////////////////////////////////////////////////////////////
		
		switch(opts.getCharVal("outMode")) {
		case 's' :
		case 'S' : printMode = PRINT_SMILES; break;
		case 'a' :
		case 'A' : printMode = PRINT_GRAPHS; break;
		case 'g' :
		case 'G' : printMode = PRINT_GML; break;
		case 'r' :
		case 'R' : printMode = PRINT_REACTIONS; break;
		case 'n' :
		case 'N' : printMode = PRINT_REACTION_NETWORK; break;
		default:
				std::ostringstream oss;
				oss	<<"output mode outmode='" <<opts.getCharVal("outMode") <<"'"
					<<" not supported";
				throw ArgException(oss.str());
		}
		
		
		//////////////////////////////////////////////////////////////
		// setup reaction rate calculation
		//////////////////////////////////////////////////////////////
		
		if (printMode == PRINT_REACTIONS || printMode == PRINT_REACTION_NETWORK) {
			switch(opts.getCharVal("rate")) {
			case 'n' : 
			case 'N' : {
				  // creates a dummy rate calculator to enforce ITS generation
				rateCalc = new RRC_TState();
				break;
			}
			case 'm' : 
			case 'M' : { 
				  // creates a web server based rate calculator using MOPAC
				RRC_QuantumMechanics *rrc = new RRC_QuantumMechanics( RRC_QuantumMechanics::QM_MOPAC );
				int errorCode = rrc->testServerAvailability();
				if (errorCode != 0) {
					std::ostringstream oss;
					oss	<<"reaction rate calculation server '"
						<<rrc->getTCPhost() <<":" <<rrc->getTCPport()
						<<"' not working properly!\n\t Error code = " 
						<<errorCode <<" = " <<strerror(errorCode);
					throw ArgException(oss.str());
				}
				rateCalc = rrc;
				break;
			}
			case 'j' : 
			case 'J' : {
				  // creates a web server based rate calculator using JAGGUAR
				RRC_QuantumMechanics *rrc = new RRC_QuantumMechanics( RRC_QuantumMechanics::QM_JAGGUAR );
				int errorCode = rrc->testServerAvailability();
				if (errorCode != 0) {
					std::ostringstream oss;
					oss	<<"reaction rate calculation server ' " 
						<<rrc->getTCPhost() <<":" <<rrc->getTCPport()
						<<"' not working properly!\n\t Error code = " 
						<<errorCode <<" = " <<strerror(errorCode);
					throw ArgException(oss.str());
				}
				rateCalc = rrc;
				break;
			}
			case 'a' :
			case 'A' : {
				  // creates Arrhenius law based rate calculator using the
				  // energy terms estimated by the molecule decomposition
				  // approach by Jankowski et al.
				ec_molDec = new EC_MoleculeDecomposition();
				  // parse Arrhenius parameter
				double arrh_kT = opts.getDoubleVal("rateArrT");
				if( arrh_kT == 0.0 ) {
					throw ArgException("Arrhenius reaction rate temperature parameter 'rateArrT' has to be non-zero");
				}
				double arrh_A = 1.0;
				rateCalc = new RRC_ArrheniusLaw( *ec_molDec, arrh_kT, arrh_A );
				break;
			}
			default  : 
				std::ostringstream oss;
				oss	<<"reaction rate calculation mode rate='" 
					<<opts.getCharVal("rate") <<"'"
					<<" not supported";
				throw ArgException(oss.str());
			}
		}
		

		//////////////////////////////////////////////////////////////
		// setup aromaticity perception model
		//////////////////////////////////////////////////////////////

		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # LOAD AROMATICITY MODELS ..."; out->flush();
		}
		switch(opts.getCharVal("aromaticity")) {
		case 'o' :
		case 'O' :
			aromaticityPrediction = new AP_NSPDK("OpenBabel:2013");
			assert( ((AP_NSPDK*)aromaticityPrediction)->getModel() != NULL );
			break;
		case 'm' :
		case 'M' :
			aromaticityPrediction = new AP_NSPDK("Marvin:general:2013");
			assert( ((AP_NSPDK*)aromaticityPrediction)->getModel() != NULL );
			break;
		case 'n' :
		case 'N' :
			aromaticityPrediction = new AP_disabled();
			break;
		default:
				std::ostringstream oss;
				oss	<<"aromaticity perception type '" <<opts.getCharVal("aromaticity") <<"'"
					<<" is not supported";
				throw ArgException(oss.str());
		}
		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}


		
		//////////////////////////////////////////////////////////////
		// parse groups from input
		//////////////////////////////////////////////////////////////


		if (opts.argExist("groups")) {
			if (infoMode == OUT_VERBOSE) {
				(*out) <<"\n # PARSE GROUPS ..."; out->flush();
			}
			  // check if input specifications are compatible
			if ( opts.argExist("smiles") && opts.getStrVal("groups") == opts.getStrVal("smiles") ) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and smiles from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}
			if ( opts.argExist("mols") && opts.getStrVal("groups") == opts.getStrVal("mols") ) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and mols from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}
			if ( opts.argExist("rules") && opts.getStrVal("groups") == opts.getStrVal("rules") ) {
				std::ostringstream oss;
				oss	<<"Cannot read both groups and rules from '" <<opts.getStrVal("groups") <<"'";
				throw ArgException(oss.str());
			}

			  // parse molecule group definitions that can be abbreviated in
			  // rules and molecules
			parseGroups( opts.getStrVal("groups"), moleculeGroups );

			if (infoMode == OUT_VERBOSE) {
				(*out) <<" DONE\n"; out->flush();
			}
		}

		//////////////////////////////////////////////////////////////
		// parse Rules from input
		//////////////////////////////////////////////////////////////
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # PARSE RULES ..."; out->flush();
		}
		std::vector<ggl::chem::ChemRule> rules;
		  // parse all rules for given input
		parseRules( opts.getStrVal("rules"), rules, moleculeGroups );

		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}
		
		if (rules.empty()) {
			(*out)	<<"\n PROBLEM : no rules found in given input!\n"
						<<std::endl;
			return 0;
		}
		
		if (!opts.argExist("noRuleCheck")) {
			if (infoMode == OUT_VERBOSE) {
				(*out) <<"\n # CHECK RULES ..."; out->flush();
			}
			  // check if all rules are valid
			bool allRulesOK = true;
			for (size_t r=0; r<rules.size(); ++r) {
				size_t conStatus = rules.at(r).isConsistent();
				if( conStatus != ggl::chem::ChemRule::C_Consistent ) {
					allRulesOK = false;
					(*out) <<"\n PROBLEM : rule " <<(r+1) <<" '"
							<<rules.at(r).getID() <<"' is not chemically correct"
								" or contains unsupported properties:\n";
					rules.at(r).decodeConsistencyStatus( conStatus, (*out) );
				}
			}
			if (!allRulesOK) {
				return 0;
			}
			if (infoMode == OUT_VERBOSE) {
				(*out) <<" DONE\n"; out->flush();
			}
		}
		
		
		  // generate left side pattern of each rule needed for its application
		  // and store separate rule lists based on the component number of
		  // their left side pattern
		for (size_t r=0; r<rules.size(); ++r) {
			ggl::chem::LeftSidePattern* pattern = new ggl::chem::LeftSidePattern(rules[r]);
			  // store in the pattern list according to the component number
			size_t compNumber = pattern->getFirstOfEachComponent().size();
			rulePattern[compNumber].push_back( pattern );
		}
		
		//////////////////////////////////////////////////////////////
		// parse molecules from input
		//////////////////////////////////////////////////////////////
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<"\n # PARSE INPUT ..."; out->flush();
		}
		if (opts.argExist("smiles")) {
			parseSMILES( opts.getStrVal("smiles"), targetSmiles, moleculeGroups );
		}
		if (opts.argExist("mols")) {
			parseMolGML( opts.getStrVal("mols"), targetSmiles, moleculeGroups );
		}
		if (infoMode == OUT_VERBOSE) {
			(*out) <<" DONE\n"; out->flush();
		}
		
		if (!opts.argExist("noInputCorrection")) {
			if (infoMode == OUT_VERBOSE) {
				(*out) <<"\n # CORRECT INPUT ..."; out->flush();
			}
			  // correct aromaticity and adjacent protons of input molecules
			correctInputMolecules( targetSmiles, producedSmiles, aromaticityPrediction, true );
			  // clear input molecules
			for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
				delete it->second;
			}
			targetSmiles.clear();
			  // setup target molecules
			std::swap( producedSmiles, targetSmiles );
			if (infoMode == OUT_VERBOSE) {
				(*out) <<" DONE\n"; out->flush();
			}
		}

		if (!opts.argExist("noInputCheck")) {
			if (infoMode == OUT_VERBOSE) {
				(*out) <<"\n # CHECK INPUT ..."; out->flush();
			}
			  // check if all rules are valid
			bool allMolOK = true;
			for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
				size_t conStatus = ggl::chem::MoleculeUtil::isConsistent( *(it->second) );
				if( conStatus != ggl::chem::MoleculeUtil::C_Consistent ) {
					allMolOK = false;
					(*out) <<"\n PROBLEM : molecule '"
							<<it->first <<"' is not chemically correct"
								" or contains unsupported properties:\n";
					ggl::chem::MoleculeUtil::decodeConsistencyStatus( conStatus, (*out) );
//					(*out) <<ggl::chem::Molecule_Graph( *(it->second) ) <<std::endl;
				}
			}
			if (!allMolOK) {
				return 0;
			}
			if (infoMode == OUT_VERBOSE) {
				(*out) <<" DONE\n"; out->flush();
			}
		}

		//////////////////////////////////////////////////////////////
		// print to stream if in verbose mode
		//////////////////////////////////////////////////////////////

		if (infoMode == OUT_VERBOSE) {

			  // print groups
			(*out) <<"\n ######## PARSED AND AVAILABLE GROUPS #######\n";
			printGroups( *out, moleculeGroups );

			(*out) <<"\n ######## PARSED RULES #######\n";
			printRules( *out, rules );
			
			for (size_t i=0; i<rules.size(); ++i ) {
			(*out) <<" transition state of "<<rules.at(i).getID();
			out->flush();
			(*out) <<" = " <<rules.at(i).getTransitionState().getSMILES() <<std::endl;
			}

			(*out) <<"\n ###### PARSED MOLECULES #####\n\n";
			printSMILES( *out, targetSmiles );
			(*out) <<std::endl; out->flush();
		}
		
		//////////////////////////////////////////////////////////////
		// perform rule application iterations
		//////////////////////////////////////////////////////////////


		ggl::Graph_Storage* gs = NULL;
		
		  // set up graph matcher
		sgm::SGM_vf2 sgm;
		
		  // the produced reactions
		MR_Reactions::Reaction_Container producedReactions;
		
		  // perform iterations
		for (size_t it=0; it<iterations; ++it ) {
			  // check if something to do
			if (targetSmiles.size() == 0)
				break;
			
			if (printMode == PRINT_REACTIONS || printMode == PRINT_REACTION_NETWORK) {
				
				(*out)	<<"\n " <<(it) <<". iteration done : molecules = "
							<<(targetSmiles.size()+producedSmiles.size()) <<std::endl;
				
				if (it==0) {// in first iteration are targetSmiles the produced ones of 'last' iteration
					std::swap(producedSmiles, targetSmiles);
				}
				  // print all new molecules
				if (opts.argExist("v") || opts.argExist("showNew")) {
					(*out) <<"\n new molecules :\n";
					  // container for sorted output
					std::vector<const std::string*> toPrint;
					  // print SMILES
					for (SMILES_container::iterator it=producedSmiles.begin(); it!=producedSmiles.end() ; ++it ) {
						toPrint.push_back(&(it->first));
					}
					  // sort SMILES to get a unique output order (needed for tests and useful anyway)
					std::sort(toPrint.begin(), toPrint.end(),compareStringPointer);
					for (std::vector<const std::string*>::const_iterator s=toPrint.begin(); s!=toPrint.end(); ++s) {
						(*out) <<" " <<**s <<"\n";
					}
					(*out) <<std::endl;
				}
				
				SMILES_container toFill;
				applyRules(rulePattern, targetSmiles, producedSmiles, toFill, producedReactions, rateCalc, allowAllIntra, *aromaticityPrediction);
				targetSmiles.insert(producedSmiles.begin(), producedSmiles.end());
				producedSmiles.clear();
				producedSmiles.insert( toFill.begin(), toFill.end() );
				toFill.clear();
				
				  // check for last iteration
				if ((it+1)==iterations) {
					targetSmiles.insert(producedSmiles.begin(), producedSmiles.end());
				}
			
			} else {
				

				(*out)	<<"\n " <<(it) <<". iteration done : molecules = "
							<<targetSmiles.size() <<std::endl;
			
				  // print all new molecules
				if (opts.argExist("v") || opts.argExist("showNew")) {
					(*out) <<"\n new molecules :\n";
					  // container for sorted output
					std::vector<const std::string*> toPrint;
					  // print SMILES
					for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
						if (producedSmiles.find(it->first) == producedSmiles.end()) {
							toPrint.push_back(&(it->first));
						}
					}
					  // sort SMILES to get a unique output order (needed for tests and useful anyway)
					std::sort(toPrint.begin(), toPrint.end(),compareStringPointer);
					for (std::vector<const std::string*>::const_iterator s=toPrint.begin(); s!=toPrint.end(); ++s) {
						(*out) <<" " <<**s <<"\n";
					}
					(*out) <<std::endl;
				}

//			// TODO does only function if all rules are applied and all old molecules are present in next iteration
//			// copy all start molecules into produced set
				for (SMILES_container::iterator oldSMILES=targetSmiles.begin(); oldSMILES!=targetSmiles.end() ; ++oldSMILES) {
					if (producedSmiles.find(oldSMILES->first)==producedSmiles.end()) {
						producedSmiles.insert(*oldSMILES);
					}
				}

				  // set up storage interface for MR_ApplyRule
				gs = new ggl::chem::GS_SMILES_MOLp<SMILES_container>(producedSmiles);
			
				  // FOR EACH set of rules with equal number of connected components
				  //          in leftsidepattern :
				  // find matches and apply rules
				for (RulePatternMap::const_iterator pat = rulePattern.begin();
						pat != rulePattern.end(); ++pat)
				{
//					size_t compNbr= pat->first;
						
					  // set up Rule applyer 
					  // check if multicomponent rule can be applied to a single molecule or not
					sgm::Match_Reporter * mr_applyRule = NULL;
	
					mr_applyRule = new ggl::MR_ApplyRule( *gs, !allowAllIntra );
//					mr_applyRule = new ggl::MR_ApplyRule( *gs, (allowAllIntra ? 1 : compNbr), false);
					
					// for all rules
					for (size_t curRule=0; curRule<pat->second.size(); ++curRule) {

						// store old target smiles size to determine the number
						// of new molecules produced by this reaction
						size_t oldProducedSmilesSize = producedSmiles.size();
						if (infoMode == OUT_VERBOSE) {
							(*out)	<<" iteration " <<(it+1)<<" :"
										<<" apply reaction '"
										<<static_cast<const ggl::chem::LeftSidePattern*>(pat->second.at(curRule))->getRule().getID()
										<<"'";
							(*out).flush();
						}

						  // set up symmetry breaking conditions for current rule
						sgm::PA_OrderCheck ga = static_cast<const ggl::chem::LeftSidePattern*>(pat->second.at(curRule))->getGraphAutomorphism();
						  // rule application
						singleRuleApplication(  sgm
												, *(static_cast<const ggl::chem::LeftSidePattern*>(pat->second.at(curRule)))
												, targetSmiles
												, *mr_applyRule
												, &ga );

						// print number of molecules produced by this reaction
						if (infoMode == OUT_VERBOSE) {
							(*out) <<" : "
										<<(producedSmiles.size()-oldProducedSmilesSize)
										<<std::endl;
						}
					}
					
					delete mr_applyRule; mr_applyRule = NULL;
				}
			
			  // make all produced SMILES target molecules for the next iteration
//			for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
//				if (producedSmiles.find(it->first) == producedSmiles.end()) {
//					producedSmiles[it->first] = it->second;
//				} else {
//					delete it->second;
//					it->second = NULL;
//				}
//			}
			
				std::swap(producedSmiles, targetSmiles);
				
				if (gs != NULL) {
					delete gs; gs = NULL;
				}
			}
			
		} // end rule application iteration loop
		
		//////////////////////////////////////////////////////////////
		// write output to stream
		//////////////////////////////////////////////////////////////
		
		(*out)	<<"\n " <<(iterations) <<". iteration done : molecules = "
					<<targetSmiles.size() <<"\n" 
					<<std::endl;
		
		  // print all new molecules
		if (opts.argExist("v") || opts.argExist("showNew")) {
			(*out) <<" new molecules :\n";
			  // container for sorted output
			std::vector<const std::string*> toPrint;
			  // print SMILES
			if (printMode == PRINT_REACTIONS || printMode == PRINT_REACTION_NETWORK) {
				for (SMILES_container::iterator it=producedSmiles.begin(); it!=producedSmiles.end() ; ++it ) {
					toPrint.push_back(&(it->first));
				}
			} else {
				for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
					if (producedSmiles.find(it->first) == producedSmiles.end()) {
						toPrint.push_back(&(it->first));
					}
				}
			}
			  // sort SMILES to get a unique output order (needed for tests and useful anyway)
			std::sort(toPrint.begin(), toPrint.end(), compareStringPointer);
			for (std::vector<const std::string*>::const_iterator s=toPrint.begin(); s!=toPrint.end(); ++s) {
				(*out) <<" " <<**s <<"\n";
			}
			(*out) <<std::endl;
		}
		
		if (infoMode == OUT_VERBOSE) {
			(*out) <<std::endl;
			
			(*out) <<"\n\n ######## FINAL SMILES ######\n\n";
		}
		  // print according output
		switch (printMode) {
			case PRINT_SMILES:
				printSMILES( *out, targetSmiles );
				break;
			case PRINT_GRAPHS:
				printMoleculeGraphs( *out, targetSmiles );
				break;
			case PRINT_GML:
				printMoleculeGraphsGML( *out, targetSmiles, true );
				break;
			case PRINT_REACTIONS:
				(*out) <<"\n number of produced reactions = " <<producedReactions.size() <<std::endl;
				(*out) <<"\n Produced Reactions :\n\n";
				for ( MR_Reactions::Reaction_Container::const_iterator 
						curReaction = producedReactions.begin();
						curReaction != producedReactions.end(); ++curReaction)
				{
					(*out) <<(*curReaction) <<"\n";
				}
				break;
			case PRINT_REACTION_NETWORK:
			{
				  // print reactions and store network information
				std::map< std::string, size_t > reaction2nodeid;
				size_t reactionID = 0;
				std::map< std::string, size_t > reactionMolecules;
				size_t moleculeID = 0;
				(*out) <<"\n number of produced reactions = " <<producedReactions.size() <<std::endl;
				(*out) <<"\n Produced Reactions :\n\n";
				for ( MR_Reactions::Reaction_Container::const_iterator
						curReaction = producedReactions.begin();
						curReaction != producedReactions.end(); ++curReaction)
				{
					(*out) <<(*curReaction) <<"\n";
					  // set new reaction ID if the reaction type is unknown so far
					if (reaction2nodeid.find(curReaction->rule_id) == reaction2nodeid.end()) {
						reaction2nodeid[curReaction->rule_id] = reactionID;
						++reactionID;
					}
					  // store all occurring molecules, ie. node in the network
					for (Reaction::Metabolite_Container::const_iterator m = curReaction->metabolites.begin();
							m != curReaction->metabolites.end(); ++m)
					{
						if (reactionMolecules.find(*m) == reactionMolecules.end()) {
							reactionMolecules[*m] = moleculeID;
							++moleculeID;
						}
					}
					for (Reaction::Product_Container::const_iterator p = curReaction->products.begin();
							p != curReaction->products.end(); ++p)
					{
						if (reactionMolecules.find(*p) == reactionMolecules.end()) {
							reactionMolecules[*p] = moleculeID;
							++moleculeID;
						}
					}
				}

				  // print reaction network in DOT format
				(*out) <<"\n\ndigraph reactionNetwork {\n";
				  // print all participating molecules
				for (std::map< std::string, size_t >::const_iterator m = reactionMolecules.begin();
						m != reactionMolecules.end(); ++m)
				{
					(*out) <<"  M" <<m->second <<" [shape=oval label=\""<<m->first<<"\"];\n";
				}
				  // print all reactions
				size_t reactionCount = 0;
				for ( MR_Reactions::Reaction_Container::const_iterator
						curReaction = producedReactions.begin();
						curReaction != producedReactions.end(); ++curReaction)
				{
					  // create reaction node
					(*out) <<"  " <<"R"<<reactionCount
							<<" [shape=box label=\"R"<<reaction2nodeid.find(curReaction->rule_id)->second;
					if (curReaction->rate == curReaction->rate) {
						(*out) <<" " <<curReaction->rate;
					}
					(*out) <<"\"]; // " <<curReaction->rule_id <<";\n";
					  // make connections
					for (Reaction::Metabolite_Container::const_iterator m = curReaction->metabolites.begin();
							m != curReaction->metabolites.end(); ++m)
					{
						(*out) <<"  M" <<reactionMolecules.find(*m)->second <<" -> R" <<reactionCount<<";\n";
					}
					for (Reaction::Product_Container::const_iterator p = curReaction->products.begin();
							p != curReaction->products.end(); ++p)
					{
						(*out) <<"  " <<"R" <<reactionCount<<" -> M" <<reactionMolecules.find(*p)->second <<";\n";
					}
					  // increase reaction node counter
					++reactionCount;
				}
				(*out) <<"} // visualization: eg. 'dot -Tpng -O GRAPHOUTPUTFILE'\n";
			}	break;
			default:
				break;
		}
		(*out)	<<std::endl;
		
		
	
	} catch (std::exception& ex) {
		(*out) <<"\n\n ERROR : " <<ex.what() <<"\n"<<std::endl;
		exitValue = -1;
	}
	
	
	//////////////////////////////////////////////////////////////
	// garbage collection
	//////////////////////////////////////////////////////////////

	if (aromaticityPrediction != NULL)	delete aromaticityPrediction;
	if (rateCalc != NULL)			delete rateCalc;
	if (ec_molDec != NULL)			delete ec_molDec;
	out = &std::cout;
	if (outFile != NULL)			delete outFile;
	for (RulePatternMap::iterator pat = rulePattern.begin();
			pat != rulePattern.end(); ++pat)
	{
		for (size_t p=0; p<pat->second.size(); ++p)
			delete pat->second[p];
	}
	for (SMILES_container::iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {
		delete it->second;
	}
	return exitValue;
}

//////////////////////////////////////////////////////////////////////////

/*!
 * Initialises allowed parameters and their default values and prepares an
 * information string for the tool.
 */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )  
{
	infoText = "\n"
		"Reads a list of molecules and chemical rules and applies the"
		"rules to generate new molecules.\n"
		"\n"
		"Rules have to be in GML format (use '-ruleExample' for an example)."
		"\n"
		"It is possible to specify the molecules in GML format as well."
		;
	
	allowedArgs.push_back(biu::COption(	
							"rules", true, biu::COption::STRING, 
							"Chemical rules input : 'STDIN' when to read from standard input, or a ':'-separated list of file names (use -ruleExample for a rule sample)"));
	allowedArgs.push_back(biu::COption(	
							"smiles", true, biu::COption::STRING, 
							"SMILES input : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(	
							"mols", true, biu::COption::STRING, 
							"Molecules input in GML format : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(	
							"groups", true, biu::COption::STRING,
							"Predefined molecule groups abbreviated in rules and molecules in GML format : 'STDIN' when to read from standard input, or a ':'-separated list of file names."));
	allowedArgs.push_back(biu::COption(
							"iter", true, biu::COption::INT, 
							"Number of rule application iterations",
							"1"));
	allowedArgs.push_back(biu::COption(	
							"allowAllIntra", true, biu::COption::BOOL, 
							"If present, all intra-molecular reactions are allowed, i.e. the application of rules with 2 or more unconnected components in the left side patter can applied to one molecule, otherwise NOT."));
	allowedArgs.push_back(biu::COption(	
							"showNew", true, biu::COption::BOOL, 
							"If present, all new molecules per iteration are printed."));
	allowedArgs.push_back(biu::COption(	
							"out", true, biu::COption::STRING, 
							"Output file name or 'STDOUT' when to write to standard output",
							"STDOUT"));
	allowedArgs.push_back(biu::COption(	
							"outMode", true, biu::COption::CHAR, 
							"Output mode : (S)MILES string, (A)djacency list, (G)ML graph representations, (R)eactions, or reaction (N)etwork",
							"S"));
	allowedArgs.push_back(biu::COption(	
							"rate", true, biu::COption::CHAR, 
							"Reaction rate calculation mode : (N)o rate calculation but transition state shown, (A)rrhenius law-based [see 'rateArrT'], (M)OPAC-based, or (J)AGGUAR-based",
							"N"));
	allowedArgs.push_back(biu::COption(	
							"rateArrT", true, biu::COption::DOUBLE,
							"Arrhenius rate generalized temperature parameter T for rate = exp(-deltaE/T); The energy difference 'deltaE' computation uses the decomposition approach by Jankowski et al. (2008)",
							"50"));
	allowedArgs.push_back(biu::COption(
							"aromaticity", true, biu::COption::CHAR,
							"The aromaticity perception model to be used : (M)arvin general model, (O)penBabel model, or (N)o aromaticity perception.",
							"M"));
	allowedArgs.push_back(biu::COption(
							"noInputCorrection", true, biu::COption::BOOL,
							"Dont correct the input molecules (aromaticity perception, proton filling, ...)"));
	allowedArgs.push_back(biu::COption(
							"noInputCheck", true, biu::COption::BOOL,
							"Dont check the input molecules for consistency (atom/bond label, ...)"));
	allowedArgs.push_back(biu::COption(
							"noRuleCheck", true, biu::COption::BOOL, 
							"Dont check the rules for consistency"));
	allowedArgs.push_back(biu::COption(	
							"help", true, biu::COption::BOOL, 
							"Displays help on all parameters"));
	allowedArgs.push_back(biu::COption(	
							"ruleExample", true, biu::COption::BOOL, 
							"Displays an example for the chemical reaction GML encoding"));
	allowedArgs.push_back(biu::COption(	
							"graphExample", true, biu::COption::BOOL, 
							"Displays an example for the molecule graph GML encoding"));
	allowedArgs.push_back(biu::COption(	
							"v", true, biu::COption::BOOL, 
							"Verbose output"));
	allowedArgs.push_back(biu::COption(	
							"version", true, biu::COption::BOOL, 
							"Version information"));
}

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraph( std::ostream& out, ggl::chem::Molecule & m )
{
	sgm::Graph_boost<ggl::chem::Molecule> g(m);
	for (size_t i=0; i<g.getNodeNumber(); ++i ) {
		out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
		for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
				eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
		{
			out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
		}
		out <<" |\n";
	}

}

//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphs( std::ostream& out, SMILES_container & smiles )
{
	for (	SMILES_container::const_iterator s = smiles.begin(); 
			s!= smiles.end(); ++s )
	{
		  // parse SMILES to graph
		std::pair<ggl::chem::Molecule,int> result 
			= ggl::chem::SMILESparser::parseSMILES(s->first);
		  // check parsing result
		if (result.second != -1) {
			std::ostringstream oss;
			oss	<<"printMoleculeGraphs : parsing error in SMILES string '"
				<<s->first
				<<"' at position " <<result.second;
			throw ArgException(oss.str());
		}
		printMoleculeGraph(out, result.first);
		out <<"\n";
	}

}


//////////////////////////////////////////////////////////////////////////

void
printMoleculeGraphsGML( std::ostream& out, SMILES_container & smiles, const bool withSpaces  )
{
	size_t i=0;
	for (	SMILES_container::const_iterator s = smiles.begin(); 
			s!= smiles.end(); ++s )
	{
		  // parse SMILES to graph
		std::pair<ggl::chem::Molecule,int> result 
			= ggl::chem::SMILESparser::parseSMILES(s->first);
		  // check parsing result
		if (result.second != -1) {
			std::ostringstream oss;
			oss	<<"printMoleculeGraphs : parsing error in SMILES string "
				<<s->first
				<<" at position " <<result.second;
			throw ArgException(oss.str());
		}
		out <<((i==0)?"":"\n")<<"# result graph " <<i <<"\n";
		ggl::Graph_GML_writer::write(out, result.first, withSpaces);
		i++;
	}
	out <<std::endl;

}


//////////////////////////////////////////////////////////////////////////

void
giveRuleGMLExample()
{
	std::cout <<"\n\
====== RULE TO ENCODE =======================\n\
\n\
  Diels-Alder reaction:\n\
\n\
  0(C)   5(C) = 4(C)               0(C) - 5(C) - 4(C) \n\
                                                      \n\
   ||                     ==>       |             |   \n\
                                                      \n\
  1(C) - 2(C) = 3(C)               1(C) = 2(C) - 3(C) \n\
\n\
======= RULE IN GML =========================\n\
\n\
rule [\n\
 ruleID \"Diels-Alder reaction\"\n\
 context [\n\
   node [ id 0 label \"C\" ]\n\
   node [ id 1 label \"C\" ]\n\
   node [ id 2 label \"C\" ]\n\
   node [ id 3 label \"C\" ]\n\
   node [ id 4 label \"C\" ]\n\
   node [ id 5 label \"C\" ]\n\
 ]\n\
 left [\n\
   edge [ source 0 target 1 label \"=\" ]\n\
   edge [ source 1 target 2 label \"-\" ]\n\
   edge [ source 2 target 3 label \"=\" ]\n\
   edge [ source 4 target 5 label \"=\" ]\n\
   constrainNoEdge [ source 0 target 4 ]\n\
   constrainNoEdge [ source 3 target 5 ]\n\
 ]\n\
 right [\n\
   edge [ source 0 target 1 label \"-\" ]\n\
   edge [ source 0 target 5 label \"-\" ]\n\
   edge [ source 1 target 2 label \"=\" ]\n\
   edge [ source 2 target 3 label \"-\" ]\n\
   edge [ source 3 target 4 label \"-\" ]\n\
   edge [ source 4 target 5 label \"-\" ]\n\
 ]\n\
]\n\
\n\
=============================================\n"
                <<std::endl;
}


//////////////////////////////////////////////////////////////////////////

void
giveGraphGMLExample()
{
	std::cout <<"\n\
====== MOLECULE TO ENCODE =======================\n\
\n\
  SMILES :  'C=CC(C)=C' \n\
\n\
======= MOLECULE IN GML =========================\n\
\n\
# molecule C=CC(C)=C \n\
graph [\n\
  node [ id 0 label \"C\" ]\n\
  node [ id 1 label \"C\" ]\n\
  node [ id 2 label \"C\" ]\n\
  node [ id 3 label \"C\" ]\n\
  node [ id 4 label \"C\" ]\n\
  edge [ source 1 target 0 label \"=\" ]\n\
  edge [ source 2 target 1 label \"-\" ]\n\
  edge [ source 3 target 2 label \"-\" ]\n\
  edge [ source 4 target 2 label \"=\" ]\n\
]\n\
\n\
==============================================\n"
                <<std::endl;
}

//////////////////////////////////////////////////////////////////////////
