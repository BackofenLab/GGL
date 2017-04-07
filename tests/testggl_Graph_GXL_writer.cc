
#include <iostream>

#include "ggl/Graph.hh"
#include "ggl/Graph_GXL_writer.hh"

#include "dataTargetGraph_1.icc"
#include "dataTargetGraph_2.icc"
#include "dataTargetGraph_3.icc"

int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"           ggl::Graph_GXL_writer  \n"
				<<"==============================================\n" 
				<<std::endl;

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_1(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		std::cout <<"\n-->Graph_GXL_writer::write(graph) : \n" <<std::endl;
		ggl::Graph_GXL_writer::write(std::cout, graph);
		std::cout <<std::endl;
	}

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_2(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		std::cout <<"\n-->Graph_GXL_writer::write(graph) : \n" <<std::endl;
		ggl::Graph_GXL_writer::write(std::cout, graph);
		std::cout <<std::endl;
	}

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_3(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		std::cout <<"\n-->Graph_GXL_writer::write(graph) : \n" <<std::endl;
		ggl::Graph_GXL_writer::write(std::cout, graph);
		std::cout <<std::endl;
	}

	std::cout	<<"\n"
				<<"==============  END TEST  ==================\n"
				<<std::endl;
	return 0;
}
