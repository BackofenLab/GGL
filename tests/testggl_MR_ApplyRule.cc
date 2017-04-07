
#include <iostream>

#include "ggl/Rule.hh"
#include "ggl/RuleGraph.hh"
#include "ggl/GS_stream.hh"
#include "ggl/MR_ApplyRule.hh"

#include "sgm/Pattern.hh"
#include "sgm/MC_Node.hh"
#include "sgm/MC_Edge.hh"
#include "sgm/MR_stream.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"
#include "sgm/Graph_boostV_p.hh"

#include "dataRule_1.icc"
#include "dataRule_2.icc"
#include "dataRule_3.icc"
#include "dataTargetGraph_1.icc"
#include "dataTargetGraph_2.icc"
#include "dataTargetGraph_3.icc"
#include "utilPrintRule.icc"

void performTest( const ggl::Rule& r
					, const int targetNumber
//					, const size_t minMatchedTargets
					, const bool addEachComponent ) {

	////////////////////////////////////////////////////////////////////////////
	
//	std::cout <<"\n=============  RULE GRAPH GRAPH  ================\n" <<std::endl;
//
//	std::cout <<"-->left side :\n" <<ggl::LeftSidePattern(r) <<std::endl;
//	std::cout <<"-->right side :\n" <<ggl::RightSidePattern(r) <<std::endl;

	std::cout <<"\n============  CREATE TARGET GRAPH  ==============\n" <<std::endl;

	if (targetNumber>0)
		std::cout <<"-->create Target = getTarget_"<<targetNumber<<"()" <<std::endl;
	else
		std::cout <<"-->create Target = getTarget_*()" <<std::endl;

	std::string graphString = "";
	ggl::Graph t;
	std::vector <const ggl::Graph*> tVec;
	
	sgm::Graph_Interface* tg = NULL;
	
	switch (targetNumber) {
	case 0: {
		tVec.push_back( new ggl::Graph(getTarget_1(graphString)));
		std::string tmp = "";
		tVec.push_back( new ggl::Graph(getTarget_2(tmp)));
		graphString += tmp;
		tg = new sgm::Graph_boostV_p<ggl::Graph>(tVec);
		break; }
	case 1: { 
		t = getTarget_1(graphString);
		tg = new sgm::Graph_boost<ggl::Graph>(t);
		break; }
	case 2: {
		t = getTarget_2(graphString);
		tg = new sgm::Graph_boost<ggl::Graph>(t);
		break; }
	case 3: {
		t = getTarget_3(graphString);
		tg = new sgm::Graph_boost<ggl::Graph>(t);
		break; }
	default:
		assert(false /* wrong target number given */);
	}
	
	std::cout <<"\n" <<graphString <<std::endl;
	
	  // usefull to print Graph objects to stream :)
	assert(tg != NULL);
	print(*tg);
	
	std::cout <<"\n================= MATCHING ======================\n" <<std::endl;
	
	
	if (r.getUsedWildcard() != NULL) {
		std::cout <<" using WILDCARD = '" <<*(r.getUsedWildcard()) <<"'\n" <<std::endl;
	} else {
		std::cout <<" using NO WILDCARD\n" <<std::endl;
	}
	ggl::LeftSidePattern ls(r);

	std::cout <<" matches :\n\n";
	
	  // check for matches
	sgm::MR_stream mr(std::cout);
	sgm::SGM_vf2 sgm;
	
	sgm.findMatches( ls, *tg, mr, UINT_MAX );

	std::cout <<"\n=================================================\n" <<std::endl;
	
	std::cout <<" results :\n\n";
	
	  // apply rule and print result to stream
	ggl::GS_stream gs(std::cout);
	ggl::MR_ApplyRule mr_rule(gs, addEachComponent);
//	ggl::MR_ApplyRule mr_rule(gs, minMatchedTargets, addEachComponent);
	
	sgm.findMatches( ls, *tg, mr_rule, UINT_MAX );

	std::cout <<"\n=================================================\n" <<std::endl;

	  // clear memory
	delete tg;
	  // clear graph container
	for (size_t i=0;i<tVec.size(); ++i)
		delete tVec[i];
	tVec.clear();
}


void performTestP( const ggl::Rule& r
					, const int targetNumber
//					, const size_t minMatchedTargets
					, const bool addEachComponent ) {

	////////////////////////////////////////////////////////////////////////////
	
	std::cout <<"\n============  CREATE TARGET GRAPH  ==============\n" <<std::endl;

	if (targetNumber>0)
		std::cout <<"-->create Target = getTarget_"<<targetNumber<<"()" <<std::endl;
	else
		std::cout <<"-->create Target = getTarget_*()" <<std::endl;

	std::string graphString = "";
	ggl::Graph t;
	std::vector < const ggl::Graph*> tVec;
	
	sgm::Graph_Interface* tg = NULL;
	
	switch (targetNumber) {
	case 0: {
		tVec.push_back( new ggl::Graph(getTarget_1(graphString)));
		std::string tmp = "";
		tVec.push_back( new ggl::Graph(getTarget_2(tmp)));
		graphString += tmp;
		tg = new sgm::Graph_boostV_p<ggl::Graph>(tVec);
		break; }
	default:
		assert(false /* wrong target number given */);
	}
	
	std::cout <<"\n" <<graphString <<std::endl;
	
	  // usefull to print Graph objects to stream :)
	assert(tg != NULL);
	print(*tg);
	
	
	std::cout <<"\n================= MATCHING ======================\n" <<std::endl;
	
	ggl::LeftSidePattern ls(r);
	
	std::cout <<" matches :\n\n";
	
	  // check for matches
	sgm::MR_stream mr(std::cout);
	sgm::SGM_vf2 sgm;
	
	sgm.findMatches( ls, *tg, mr, UINT_MAX );

	std::cout <<"\n=================================================\n" <<std::endl;
	
	std::cout	
//				<<" minMatchedTargets : " <<minMatchedTargets <<"\n"
				<<" addEachComponent  : " <<(addEachComponent?"true":"false") <<"\n"
				<<std::endl;
	
	std::cout <<" results :\n\n";
	
	  // apply rule and print result to stream
	ggl::GS_stream gs(std::cout);
	ggl::MR_ApplyRule mr_rule(gs, addEachComponent);
//	ggl::MR_ApplyRule mr_rule(gs, minMatchedTargets, addEachComponent);
	
	sgm.findMatches( ls, *tg, mr_rule, UINT_MAX );

	std::cout <<"\n=================================================\n" <<std::endl;

	  // clear memory
	delete tg;
	if (targetNumber == 0) {
		for (size_t i=0; i<tVec.size(); i++) {
			delete tVec[i];
		}
	}

}

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             ggl::MR_ApplyRule  \n" 
				<<"==============================================\n" 
				<<std::endl;
	

	std::string ruleString;
	std::string wildcard = "*";
	
	  // run rule 1 
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;
		
		std::cout <<"-->create Rule = getRule_1()\n" <<std::endl;
	
		ggl::Rule r(getRule_1(ruleString),"rule1");
		
		std::cout <<"\n" <<ruleString <<std::endl;
	
		  // apply rule to Target
		performTest( r, 1, false);
		performTest( r, 2, false );
	}

	  // run rule 1 WILDCARD
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;

		std::cout <<"-->create Rule = getRule_1()" <<std::endl;

		ggl::Rule r(getRule_1_WC(ruleString),"rule1-wc");

		std::cout <<"\n" <<ruleString <<std::endl;

		  // apply rule to Target
		performTest( r, 1, false );

		r.setUsedWildcard(&wildcard);
		performTest( r, 1, false );
		r.setUsedWildcard(NULL);
	}

	  // run rule 1 with constraints
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;

		std::cout <<"-->create RuleConstraints:\n" <<std::endl;

		std::vector< sgm::Pattern_Interface::Match_Constraint* > rC;
		rC.push_back( new sgm::MC_NodeAdjacency( 2, sgm::MC_NodeAdjacency::MC_EQ, 0, "B") );
		std::cout <<" RC( node(2), =, 0, 'B' )" <<std::endl;

		std::cout <<"\n-->create Rule = getRule_1()" <<std::endl;

		ggl::Rule r(getRule_1(ruleString), "rule1", rC);

		std::cout <<"\n" <<ruleString <<std::endl;


		  // apply rule to Target
		performTest( r, 1, false);

		delete (*rC.begin());
	}

	  // run rule 2
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;
		
		std::cout <<"-->create Rule = getRule_2()" <<std::endl;
	
		ggl::Rule r(getRule_2(ruleString),"rule2");
		
		std::cout <<"\n" <<ruleString <<std::endl;
	
		  // apply rule to Target
		performTest( r, 1, false );
		performTest( r, 2, false );
	}

	  // run rule 2 with constraints
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;

		std::cout <<"-->create RuleConstraints:\n" <<std::endl;

		std::vector< sgm::Pattern_Interface::Match_Constraint * > rC;
		std::set<std::string> allowedLabel;
		rC.push_back( new sgm::MC_NodeAdjacency( 1, sgm::MC_NodeAdjacency::MC_L, 3, allowedLabel, allowedLabel) );
		std::cout <<" RC( node(1), <, 3, any, any )" <<std::endl;

		std::cout <<"\n-->create Rule = getRule_2()" <<std::endl;

		ggl::Rule r(getRule_2(ruleString), "rule2", rC);

		std::cout <<"\n" <<ruleString <<std::endl;

		  // apply rule to Target
		performTest( r, 2, false );

		delete (*rC.begin());
	}

	  // run rule 2 with constraints
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;

		std::cout <<"-->create RuleConstraints:\n" <<std::endl;

		std::vector< sgm::Pattern_Interface::Match_Constraint * > rC;
		std::set<std::string> nodeLabel;
		std::set<std::string> edgeLabel;
		edgeLabel.insert("-2-");
		rC.push_back( new sgm::MC_NodeAdjacency( 1, sgm::MC_NodeAdjacency::MC_G, 0, nodeLabel, edgeLabel) );
		std::cout <<" RC( node(2), >, 0, any, '-2-' )" <<std::endl;

		std::cout <<"\n-->create Rule = getRule_2()" <<std::endl;

		ggl::Rule r(getRule_2(ruleString), "rule2", rC);

		std::cout <<"\n" <<ruleString <<std::endl;

		  // apply rule to Target
		performTest( r, 2, false );

		delete (*rC.begin());
	}

	  // run rule 3 
	{
		std::cout <<"\n================  CREATE RULE  ==================\n" <<std::endl;
		
		std::cout <<"-->create Rule = getRule_3()" <<std::endl;
	
		ggl::Rule r(getRule_3(ruleString),"rule3");
		
		std::cout <<"\n" <<ruleString <<std::endl;
	
		  // apply rule to Target
		performTest( r, 1, false );
		performTest( r, 2, false );
		performTest( r, 3, false );
		performTest( r, 0, false );
		performTest( r, 0, true );
		performTestP( r, 0, false );
//		performTestP( r, 0, 2, false );
		performTestP( r, 0, true );
	}

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
