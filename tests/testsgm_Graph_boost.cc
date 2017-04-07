

#include <iostream>

#include "sgm/Graph_boost.hh"

#include "dataSGMgraphs_1.icc"

#include "utilPrintGraph_Interface.icc"



	void
	printIteration( const sgm::Graph_Interface& g )
	{
		for (size_t i=0; i<g.getNodeNumber(); i++) {
			std::cout <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i);
					e != g.getOutEdgesEnd(i); ++e)
			{
				std::cout <<" | " <<(*e).getToIndex() <<" (" <<(*e).getEdgeLabel() <<")";
			}
			std::cout <<" |\n";
		}

	}


int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             sgm::Graph_boost  \n" 
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
		
		std::cout <<"--> create sgm::Graph_boost \n" <<std::endl;
		GB pg( pattern );
		
		printIteration( pg );
		
		std::cout <<"\n == TARGET ==\n" <<std::endl;

		MyGraph target = getTarget_1( targetString );

		std::cout <<targetString <<std::endl;
		
		std::cout <<"--> create sgm::Graph_boost \n" <<std::endl;
		GB tg( target );
		
		printIteration( tg );
		
		std::cout <<"\n--> comparison : (pattern == pattern) = "; 
		std::cout.flush();
		std::cout <<(pg==pg?"true":"false") <<std::endl;
		std::cout <<"\n--> comparison : (pattern != pattern) = "; 
		std::cout.flush();
		std::cout <<(pg!=pg?"true":"false") <<std::endl;
		std::cout <<"\n--> comparison : (pattern == target) = "; 
		std::cout.flush();
		std::cout <<(pg==tg?"true":"false") <<std::endl;
		std::cout <<"\n--> comparison : (pattern != target) = "; 
		std::cout.flush();
		std::cout <<(pg!=tg?"true":"false") <<std::endl;
	}

	{
		std::cout <<"\n == GRAPH ==\n" <<std::endl;

		std::string graph = "graph [ \n"
				" node [ id 1 label A ] \n"
				" node [ id 2 label B ] \n"
				" edge [ source 1 target 2 label \"e1\" ] \n"
				" edge [ source 1 target 2 label \"e2\" ] \n"
				"]";
		std::cout <<" Graph in GML = " << graph <<std::endl;


		  // the target graph to fill
		MyGraph t;

		MyGraph::vertex_descriptor v;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), t );

		MyGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type
			edgeLabel = boost::get( boost::edge_name_t(), t );

		  // add vertices
		v = boost::add_vertex(t); nodeLabel[v] = "A";
		v = boost::add_vertex(t); nodeLabel[v] = "B";
		  // add edges
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "e1";
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "e2";

		std::cout <<" number of edges = " <<boost::num_edges(t)<<std::endl;

		std::cout <<"--> create sgm::Graph_boost \n" <<std::endl;
		GB tg1( t );

		printIteration( tg1 );



	}

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}

