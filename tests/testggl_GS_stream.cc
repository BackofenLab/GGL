

#include <iostream>

#include "ggl/GS_stream.hh"

#include "dataTargetGraph_1.icc"

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"              ggl::GS_stream  \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	std::string graphString = "";
	
	std::cout <<"-->create ggl::Graph = getTarget_1()" <<std::endl;
	
	ggl::Graph t = getTarget_1(graphString);

	
	std::cout <<"\n" <<"filled graph is :\n\n" <<graphString <<std::endl;
	
	std::cout <<"-->print graph via GS_stream.add()" <<std::endl;
	  // usefull to print Graph objects to stream :)
	ggl::GS_stream gs(std::cout);
	gs.add(t);

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
