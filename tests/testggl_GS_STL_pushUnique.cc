
#include <iostream>

#include "ggl/Rule.hh"
#include "ggl/GS_STL.hh"


#include "sgm/MR_stream.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/Graph_boost.hh"

#include "sgm/GM_vf2.hh"

#include "dataTargetGraph_1.icc"
#include "dataTargetGraph_2.icc"
#include "dataTargetGraph_3.icc"
#include "utilPrintRule.icc"



int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"             ggl::GS_STL_pushUnique  \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	std::cout <<"-->create std::vector<Graph> storage" <<std::endl; 
	std::vector<ggl::Graph> storage;
	std::cout <<"-->create sgm::GM_vf2 matcher" <<std::endl; 
	sgm::GM_vf2 matcher;
	std::cout <<"-->create ggl::GS_STL_pushUnique gs(storage,matcher)" <<std::endl; 
	ggl::GS_STL_pushUnique gs(storage, matcher);

	
	std::cout <<"-->create ggl::Graph = getTarget_1()" <<std::endl;
	
	std::string graphString1 = "";
	ggl::Graph t1 = getTarget_1(graphString1);
	std::string graphString2 = "";
	ggl::Graph t2 = getTarget_2(graphString2);
	std::string graphString3 = "";
	ggl::Graph t3 = getTarget_3(graphString3);
	
	std::cout <<"\n" <<" graph 1 :\n\n" <<graphString1 <<std::endl;
	std::cout <<"\n" <<" graph 2 :\n\n" <<graphString2 <<std::endl;
	std::cout <<"\n" <<" graph 3 :\n\n" <<graphString3 <<std::endl;

	std::cout <<"->gs.add(t1)" <<std::endl;
	gs.add(t1);
	std::cout <<"->gs.add(t2)" <<std::endl;
	gs.add(t2);
	std::cout <<"->gs.add(t1)" <<std::endl;
	gs.add(t1);
	std::cout <<"->gs.add(t3)" <<std::endl;
	gs.add(t3);
	std::cout <<"->gs.add(t2)" <<std::endl;
	gs.add(t2);
	std::cout <<"->gs.add(t2)" <<std::endl;
	gs.add(t2);

	std::cout <<"\n-->resulting storage :\n" <<std::endl;
	for (size_t i=0; i<storage.size(); i++) {
		sgm::Graph_boost<ggl::Graph> gr(storage[i]);
		print(gr);
		std::cout <<std::endl;
	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
