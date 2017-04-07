
#include <iostream>

#include "ggl/Graph.hh"
#include "ggl/Graph_GML_writer.hh"

#include "dataTargetGraph_1.icc"
#include "dataTargetGraph_2.icc"
#include "dataTargetGraph_3.icc"

int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"           ggl::Graph_GML_writer  \n" 
				<<"==============================================\n" 
				<<std::endl;

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_1(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		std::cout <<"\n-->Graph_GML_writer::write(graph,true) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(std::cout, graph, true);
		std::cout <<std::endl;
		std::cout <<"\n-->Graph_GML_writer::write(graph,false) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(std::cout, graph, false);
		std::cout <<std::endl;
		std::cout <<"\n-->operator<<(graph) : \n" <<std::endl;
		std::cout <<graph <<std::endl;
	}

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_2(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		std::cout <<"\n-->Graph_GML_writer::write(graph,true) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(std::cout, graph, true);
		std::cout <<std::endl;
		std::cout <<"\n-->Graph_GML_writer::write(graph,false) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(std::cout, graph, false);
		std::cout <<std::endl;
	}

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_3(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		std::cout <<"\n-->Graph_GML_writer::write(graph,true) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(std::cout, graph, true);
		std::cout <<std::endl;
		std::cout <<"\n-->Graph_GML_writer::write(graph,false) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(std::cout, graph, false);
		std::cout <<std::endl;
	}

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<std::endl;
	return 0;
}
