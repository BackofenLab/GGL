

#include <iostream>

#include "sgm/Graph_boost.hh"
#include "sgm/Graph_NodeLabelPrefix.hh"

#ifndef DATASGMGRAPHS_GRAPH_CLASS
#define DATASGMGRAPHS_GRAPH_CLASS

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

	//! The properties available for the nodes of a graph
	typedef	boost::property<	boost::vertex_index_t, size_t
			, boost::property<	boost::vertex_name_t,  std::string
				> >
					MyGraph_NodeProperties;

	//! The properties available for the edges of a graph
	typedef	boost::property<	boost::edge_name_t, std::string
				>
					MyGraph_EdgeProperties;

	//! The definition of undirected graphs
	typedef boost::adjacency_list<
					boost::vecS,      				// store edges
					boost::vecS,       				// store vertices
					boost::undirectedS,				// is an undirected graph
					MyGraph_NodeProperties,  		// (atom symbols etc)
					MyGraph_EdgeProperties   		// (edge symbols etc)
				>
					MyGraph;


#endif // DATASGMGRAPHS_GRAPH_CLASS



#include "utilPrintGraph_Interface.icc"



	MyGraph
	getNodes_1() {

		  // the target graph to fill
		MyGraph t;

		MyGraph::vertex_descriptor v;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), t );

		MyGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type
			edgeLabel = boost::get( boost::edge_name_t(), t );

		v = boost::add_vertex(t); nodeLabel[v] = "N";
		v = boost::add_vertex(t); nodeLabel[v] = "C";
		v = boost::add_vertex(t); nodeLabel[v] = "C+2";
		v = boost::add_vertex(t); nodeLabel[v] = "C:1";
		v = boost::add_vertex(t); nodeLabel[v] = "C-1:1";
		v = boost::add_vertex(t); nodeLabel[v] = "CH";
		v = boost::add_vertex(t); nodeLabel[v] = "CH:823";
		v = boost::add_vertex(t); nodeLabel[v] = "CH-1";
		v = boost::add_vertex(t); nodeLabel[v] = "CH-3:34";
		v = boost::add_vertex(t); nodeLabel[v] = "CH2";
		v = boost::add_vertex(t); nodeLabel[v] = "CH2:2";
		v = boost::add_vertex(t); nodeLabel[v] = "CH2-3";
		v = boost::add_vertex(t); nodeLabel[v] = "CH2-3:1";

		return t;
	}



int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             sgm::Graph_NodeLabelPrefix  \n"
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
		std::cout <<"\n == molecule : ignore class label ==\n" <<std::endl;

		const std::string atomNoClass = ":";
		std::cout <<"--> applied node pattern : '"<<atomNoClass<<"'\n" <<std::endl;
		
		// create graph to wrap
		MyGraph origGraph = getNodes_1();
		
		std::cout <<"--> create sgm::Graph_boost og \n" <<std::endl;
		GB og( origGraph );
		std::cout <<"--> create sgm::Graph_SubNodeLabel( og, nodePattern )\n" <<std::endl;
		sgm::Graph_NodeLabelPrefix ng( og, atomNoClass );
		
		// compare labels
		for (size_t i=0; i<og.getNodeNumber(); i++) {
			std::cout <<" #"<<i<<"\t: "<<og.getNodeLabel(i)<<"\t=> "<<ng.getNodeLabel(i)<<std::endl;
		}
		
	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}

