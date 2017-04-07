

#include <iostream>

#include "sgm/Graph_boostV_p.hh"

#include "dataSGMgraphs_1.icc"

#include "utilPrintGraph_Interface.icc"


int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             sgm::Graph_boostV_p  \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	typedef sgm::Graph_boostV_p<	MyGraph
								, boost::vertex_name_t
								, boost::edge_name_t
								, boost::vertex_index_t
							> GBV;

	std::string graphString;
	
	  // example 1
	{
		std::cout <<"\n == GRAPH ==\n" <<std::endl;
		
		std::vector< MyGraph > graphs;
		
		graphs.push_back( getPattern_1( graphString ) );
		std::cout <<graphString <<std::endl;
		graphs.push_back( getTarget_1( graphString ) );
		std::cout <<graphString <<std::endl;
		
		std::vector< const MyGraph* > graphP(graphs.size(), NULL);
		for (size_t i=0; i<graphs.size(); i++)
			graphP[i] = &(graphs[i]); 
		
		std::cout <<"--> create sgm::Graph_boostV_p g \n" <<std::endl;
		GBV g( graphP );
		
		print( g );
		
		sgm::Graph_Interface::CompLabel label;
		std::cout <<"--> connectedComponents( g ) = "; std::cout.flush();
		std::cout	<<sgm::Graph_Interface::connectedComponents( g, label ) 
					<<std::endl;
		
		std::cout <<"--> component labeling :" <<std::endl;
		for (size_t i=0; i<label.size(); i++)
			std::cout <<" " <<i <<" = " <<label[i] <<std::endl;
		

		
	}

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}

