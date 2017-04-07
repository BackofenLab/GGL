

#include <iostream>

#include "sgm/MR_Counting.hh"
#include "sgm/Graph_boost.hh"

#include "dataSGMgraphs_1.icc"

#include "utilPrintGraph_Interface.icc"



int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"              sgm::MR_Counting  \n" 
				<<"==============================================\n" 
				<<std::endl;
	

	std::cout <<"--> create sgm::MR_Counting mr\n" <<std::endl;
	sgm::MR_Counting mr;
	
	typedef sgm::Graph_boost<	MyGraph
								, boost::vertex_name_t
								, boost::edge_name_t
								, boost::vertex_index_t
							> GB;

	  // create graphs
	std::cout <<"--> create dummy graphs \n" <<std::endl;
	std::string patternString;
	std::string targetString;
						
	GB patternGraph (getPattern_1( patternString ));
	sgm::Pattern pattern( patternGraph );
	GB target (getTarget_1( targetString ));
	
	  // create dummy match:
	std::cout <<"--> create dummy match :\n" <<std::endl;
	sgm::Match match;
	match.push_back(4);
	match.push_back(3);
	match.push_back(2);
	match.push_back(5);
	for (size_t i=0; i<match.size(); i++) {
		std::cout <<" " <<match[i];
	}
	std::cout <<std::endl;
	
	std::cout <<"--> mr.reportHit(..)" <<std::endl;
	mr.reportHit(pattern, target, match);
	std::cout <<"--> mr.reportHit(..)" <<std::endl;
	mr.reportHit(pattern, target, match);
	std::cout <<"--> mr.getHits() = " <<mr.getHits() <<std::endl;
	
	std::cout <<"--> mr.resetHits(..)" <<std::endl;
	mr.resetHits();
	std::cout <<"--> mr.getHits() = " <<mr.getHits() <<std::endl;
	
	std::cout <<"--> mr.reportHit(..)" <<std::endl;
	mr.reportHit(pattern, target, match);
	std::cout <<"--> mr.getHits() = " <<mr.getHits() <<std::endl;

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}

