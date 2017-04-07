
#include <iostream>
#include <sstream>

#include "ggl/Graph.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/Graph_GML_grammar.hh"

#include "sgm/Graph_boost.hh"

#include "dataTargetGraph_1.icc"
#include "dataTargetGraph_4.icc"

#include "utilPrintGraph_Interface.icc"

int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"           ggl::Graph_GML_grammar  \n" 
				<<"==============================================\n" 
				<<std::endl;

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_1(graphString);
		
		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;
		
		std::stringstream ss;
		std::cout <<"\n-->Graph_GML_writer::write(ss,graph,true) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(ss, graph, true);
		
		std::cout <<"\n-->Graph_GML_grammar::parseGraph(ss.str())" <<std::endl;
		
		std::pair<ggl::Graph, int> ret = ggl::Graph_GML_grammar::parseGraph(ss.str());
		
		if (ret.second < 0 ) 
		{
			sgm::Graph_boost<ggl::Graph> gi(ret.first);
			std::cout <<"\n result graph : " <<std::endl;
			print(gi);
			std::cout <<std::endl;
		} else {
			std::cout <<"\n PARSING ERROR : at input position " <<ret.second 
					<<" = '" <<ss.str().at(ret.second) <<"' "
					<<" within \n'"
					<<ss.str().substr(std::max(0,((int)ret.second)-10),20)
					<<"'\n"
					<<std::endl;
		}
	}

	  // write graph
	{
		std::string graphString;
		ggl::Graph graph = getTarget_4(graphString);

		std::cout <<"\n graph to write:\n" <<graphString <<std::endl;

		std::stringstream ss;
		std::cout <<"\n-->Graph_GML_writer::write(ss,graph,true) : \n" <<std::endl;
		ggl::Graph_GML_writer::write(ss, graph, true);

		std::cout <<"\n-->Graph_GML_grammar::parseGraph(ss.str())" <<std::endl;

		std::pair<ggl::Graph, int> ret = ggl::Graph_GML_grammar::parseGraph(ss.str());

		if (ret.second < 0 )
		{
			sgm::Graph_boost<ggl::Graph> gi(ret.first);
			std::cout <<"\n result graph : " <<std::endl;
			print(gi);
			std::cout <<std::endl;
		} else {
			std::cout <<"\n PARSING ERROR : at input position " <<ret.second
					<<" = '" <<ss.str().at(ret.second) <<"' "
					<<" within \n'"
					<<ss.str().substr(std::max(0,((int)ret.second)-10),20)
					<<"'\n"
					<<std::endl;
		}
	}


	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<std::endl;
	return 0;
}
