
#include <iostream>

#include "ggl/Rule_GMLparser.hh"
#include "ggl/GS_stream.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESparser.hh"
#include "ggl/chem/SMILESwriter.hh"
#include "ggl/chem/ChemRule.hh"
#include "ggl/chem/MR_Reactions.hh"

#include "sgm/MR_stream.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"
#include "sgm/Graph_boostV_p.hh"

#include "dataChemRule_1.icc"
#include "utilPrintRule.icc"

using namespace ggl;
using namespace ggl::chem;

void performTest(	const std::vector<std::string> & inputSMILES
					, const ggl::chem::ChemRule& r
					, const bool addEachComponent
					, const ggl::chem::ReactionRateCalculation* rateCalc ) {

	////////////////////////////////////////////////////////////////////////////
	
	// create data structures needed for call
	Smiles2GraphMap knownSmiles2graph;
	std::vector <const ggl::chem::Molecule*> mVec;
	Smiles2GraphMap newSmiles2graph;
	MR_Reactions::Reaction_Container newReactions;

	std::cout <<"\n===========  CREATE MOLECULE INPUT  =============\n" <<std::endl;

	for (std::vector<std::string>::const_iterator it=inputSMILES.begin(); it!=inputSMILES.end(); it++) {
		  // parse SMILES to graph
		std::pair<ggl::chem::Molecule,int> result
			= ggl::chem::SMILESparser::parseSMILES(*it);
		  // check parsing result
		if (result.second != -1) {
			std::cout	<<"parsing error in SMILES string '"
						<<*it
						<<"' at position " <<result.second;
			return;
		}
		  // add protons
		ggl::chem::MoleculeUtil::fillProtons(result.first);
		  // derive canonical SMILES with own writer
		const std::string canSMILES =
			ggl::chem::SMILESwriter::getSMILES( result.first, true );
		  // store parsed Molecule graph
		knownSmiles2graph[canSMILES] = new ggl::chem::Molecule(result.first);
		std::cout <<" storing : '" <<canSMILES <<"'" <<std::endl;
		mVec.push_back(knownSmiles2graph[canSMILES]);
	}

	std::cout <<"\n============  CREATE TARGET GRAPH  ==============\n" <<std::endl;


	std::string graphString = "";
	sgm::Graph_boostV_p< ggl::chem::Molecule > targets(mVec);
	
	print(targets);
	
	
	std::cout <<"\n================= MATCHING ======================\n" <<std::endl;
	
	ggl::chem::LeftSidePattern ls(r);
	std::cout <<" pattern graph = \n" <<ls <<"\n" <<std::endl;
	
	std::cout <<" matches :\n\n";
	
	  // check for matches
	sgm::MR_stream mr(std::cout);
	sgm::SGM_vf2 sgm;
	
	sgm.findMatches( ls, targets, mr, UINT_MAX );

	std::cout <<"\n============== REACTION GENERATION ==============\n" <<std::endl;
	
	std::cout <<" create MR_Reactions object( addEachComponent=" <<(addEachComponent?"true":"false")<<" )" <<std::endl;
	
	  // apply rule and print result to stream
	MR_Reactions mr_reactions( knownSmiles2graph, newSmiles2graph, newReactions, addEachComponent, rateCalc);
	
	std::cout <<" search matches and create Reaction objects ..." <<std::endl;
	try {
		sgm.findMatches( ls, targets, mr_reactions, UINT_MAX );
	} catch (std::exception& ex) {
		std::cout <<"\n EXCEPTION RAISED : "<< ex.what() <<std::endl;
	}

	std::cout <<" number of reactions found : " <<newReactions.size() <<std::endl;

	std::cout <<" reactions produced :\n" <<std::endl;
	for (MR_Reactions::Reaction_Container::const_iterator it=newReactions.begin(); it!=newReactions.end(); it++) {
		std::cout <<*it <<std::endl;
	}

	std::cout <<"\n=================================================\n" <<std::endl;

	  // clear memory
	for (size_t i=0; i<mVec.size(); i++)
		delete mVec[i];

}



int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             ggl::chem::MR_Reactions  \n"
				<<"==============================================\n" 
				<<std::endl;
	

	std::string ruleString;
	


	ggl::chem::RRC_TState rateCalc_TState;
	ggl::chem::RRC_Fixed rateCalc_fixed(1.23456789);

	  // run rule 1 
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;
		
		std::pair<ggl::Rule, int> ret = ggl::Rule_GMLparser::parseRule( ruleDielsAlderGML );
		
		std::cout <<" result of rule GML parsing : " <<ret.second <<std::endl;
		
		std::cout <<" rule ID : " <<ret.first.getID() <<std::endl;
		
		ggl::chem::ChemRule chemRule(ret.first);

		std::cout <<"\n============  CREATE SMILES INPUT  ==============\n" <<std::endl;

		std::vector<std::string>  inputSMILES;
		inputSMILES.push_back( "C=CO" );
		inputSMILES.push_back( "C=CC(=C)C" );
		inputSMILES.push_back( "O=C1C=CC(=O)O1" );
		inputSMILES.push_back( "C1=CCCC=C1" );
		inputSMILES.push_back( "C1=CC=C1" );

		std::cout <<"\n input SMILES :" <<std::endl;
		for (std::vector<std::string>::const_iterator it=inputSMILES.begin(); it!=inputSMILES.end(); it++)
			std::cout <<*it <<std::endl;
		
	
		  // apply rule to Target
		performTest( inputSMILES, chemRule, false, NULL );
		performTest( inputSMILES, chemRule, false, &rateCalc_TState );
		performTest( inputSMILES, chemRule, true, &rateCalc_fixed );
	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
