
#include <iostream>

#include "sgm/Graph_boost.hh"
#include "sgm/MR_stream.hh"
#include "sgm/SGM_vf2.hh"

#include "dataSGMgraphs_1.icc"
#include "dataSGMgraphs_2.icc"
#include "dataSGMgraphs_3.icc"

#include "utilPrintGraph_Interface.icc"

	typedef sgm::Graph_boost<	MyGraph
								, boost::vertex_name_t
								, boost::edge_name_t
								, boost::vertex_index_t
							> GB;

	void
	performTest (	const MyGraph& pattern,
					const MyGraph& target,
					const std::string wildcard = std::string("") )
	{
		
		GB pattern_(pattern);
		sgm::Pattern patternGraph( pattern_ );
		if (wildcard.size() > 0) {
			patternGraph = sgm::Pattern( pattern_, wildcard );
		}

		GB targetGraph(target);

		std::cout <<"\n--> Pattern graph:\n" <<patternGraph <<std::endl;
		std::cout <<"\n--> Target graph:\n" <<targetGraph <<std::endl;

		std::cout <<"\n--> ONE pattern and ONE target\n" <<std::endl;
		
		if (wildcard.size() == 0) {
			std::cout <<"\n find ALL matches with NO_WILDCARD :\n" <<std::endl;
		} else {
			std::cout <<"\n find ALL matches with WILDCARD = '"<<wildcard<<"' :\n" <<std::endl;
		}

		std::cout <<"--> create sgm::MR_stream\n" <<std::endl;
		
		sgm::MR_stream mr(std::cout);
		size_t hits = 0;
		
		sgm::SubGraphMatching* sgm = new sgm::SGM_vf2();
		
		hits = sgm->findMatches( patternGraph, targetGraph, mr, UINT_MAX );

		std::cout <<"\n found " <<hits <<" matches\n" <<std::endl;
		

		std::cout <<"\n--> TWO pattern and ONE target\n" <<std::endl;
		
		std::vector< const sgm::Pattern_Interface* > patternGraphs;
		patternGraphs.push_back(&patternGraph);
		patternGraphs.push_back(&patternGraph);
		std::vector<sgm::Match_Reporter*> outputs;
		outputs.push_back(&mr);
		outputs.push_back(&mr);
		
		if (wildcard.size() == 0) {
			std::cout <<"\n find ALL matches with NO_WILDCARD :\n" <<std::endl;
		} else {
			std::cout <<"\n find ALL matches with WILDCARD = '"<<wildcard<<"' :\n" <<std::endl;
		}
		
		hits = sgm->findMatches( patternGraphs, targetGraph, outputs, UINT_MAX );

		std::cout <<"\n found " <<hits <<" matches\n" <<std::endl;
		
		if (wildcard.size() == 0) {
			std::cout <<"\n find ONE match with NO_WILDCARD :\n" <<std::endl;
		} else {
			std::cout <<"\n find ONE match with WILDCARD = '"<<wildcard<<"' :\n" <<std::endl;
		}

		hits = sgm->findMatches( patternGraphs, targetGraph, outputs, 1 );

		std::cout <<"\n found " <<hits <<" matches\n" <<std::endl;
		
		delete sgm;

	}

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"               sgm::SGM_vf2  \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	typedef sgm::Graph_boost<	MyGraph
								, boost::vertex_name_t
								, boost::edge_name_t
								, boost::vertex_index_t
							> GB;

	std::string patternString;
	std::string targetString;
	
	  // example 1
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_1( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_1( targetString );
		std::cout <<targetString <<std::endl;
		

		performTest( pattern, target );
	}

	  // example 2 : wildcards
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_1( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getPattern_1_WC( targetString );
		std::cout <<targetString <<std::endl;


		performTest( pattern, target );
		performTest( pattern, target, "*" );
	}

	  // example 3 : wildcards
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_1_WC( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getPattern_1( targetString );
		std::cout <<targetString <<std::endl;


		performTest( pattern, target );
		performTest( pattern, target, "*" );
	}

	  // parallel edge example 1
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_2( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_2( targetString );
		std::cout <<targetString <<std::endl;

		performTest( pattern, target );
	}


	  // example 2 : wildcards
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_2_WC( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_2( targetString );
		std::cout <<targetString <<std::endl;


		performTest( pattern, target );
		performTest( pattern, target, "*" );
	}

	  // loop edge example 1
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_3( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_3( targetString );
		std::cout <<targetString <<std::endl;

		performTest( pattern, target );
	}


	  // example 2 : wildcards
	{
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_3_WC( patternString );
		std::cout <<patternString <<std::endl;

		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_3( targetString );
		std::cout <<targetString <<std::endl;


		performTest( pattern, target );
		performTest( pattern, target, "*" );
	}



	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}


