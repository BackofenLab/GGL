

#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>

#include "sgm/SubGraph.hh"
#include "sgm/Graph_boost.hh"

#include "dataSGMgraphs_1.icc"

#include "utilPrintGraph_Interface.icc"


int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             sgm::SubGraph  \n" 
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
		//////////////////////////////////////////////////////////////////////
		
		std::cout <<"\n == PATTERN ==\n" <<std::endl;

		MyGraph pattern = getPattern_1( patternString );
		
		std::cout <<patternString <<std::endl;
		
		std::cout <<"--> create sgm::Graph_boost pg \n" <<std::endl;
		GB pg( pattern );
		
		print( pg );
		
		// set up subset of nodes
		std::cout <<"--> node subset pn : ";
		sgm::SubGraph::NodeList pn;
		pn.insert(pn.end(), 3);
		pn.insert(pn.end(), 0);
		pn.insert(pn.end(), 2);
		
		std::copy(	pn.begin(), 
					pn.end(), 
					std::ostream_iterator<size_t>(std::cout, ", "));
		std::cout <<std::endl;
		
		std::cout <<"--> create sgm::SubGraph( pg, pn)" <<std::endl;
		
		sgm::SubGraph pSub( pg, pn );
		print( pSub );
		
		//////////////////////////////////////////////////////////////////////
		
		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_1( targetString );

		std::cout <<targetString <<std::endl;
		
		std::cout <<"--> create sgm::Graph_boost tg \n" <<std::endl;
		GB tg( target );
		
		print( tg );
		
		// set up subset of nodes
		std::cout <<"--> node subset tn : ";
		sgm::SubGraph::NodeList tn;
		tn.insert(tn.end(), 3);
		tn.insert(tn.end(), 0);
		tn.insert(tn.end(), 2);
		tn.insert(tn.end(), 5);
		tn.insert(tn.end(), 7);
		tn.insert(tn.end(), 9);
		
		std::copy(	tn.begin(), 
					tn.end(), 
					std::ostream_iterator<size_t>(std::cout, ", "));
		std::cout <<std::endl;
		
		std::cout <<"--> create sgm::SubGraph( tg, tn)" <<std::endl;
		
		sgm::SubGraph tSub( tg, tn );
		print( tSub );
		
		
	}

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}

