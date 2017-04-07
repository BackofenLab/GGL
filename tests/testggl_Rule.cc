
#include <iostream>

#include "ggl/Rule.hh"

#include "dataRule_1.icc"
#include "dataRule_2.icc"
#include "dataRule_3.icc"
#include "dataRule_4.icc"
#include "utilPrintRule.icc"

void checkConnectedComponents( ggl::Rule& r)
{

	std::cout <<"\n == CHECK CONSISTENCY ==\n" <<std::endl;
	std::cout <<"-> isConsistent() = ";
	const size_t consistencyCode = r.isConsistent();
	std::cout <<(consistencyCode==ggl::Rule::C_Consistent?"true":"false") <<std::endl;
	if (consistencyCode != ggl::Rule::C_Consistent) {
		r.decodeConsistencyStatus( consistencyCode, std::cout );
	}

	sgm::Graph_Interface::CompLabel label;
	
	label.clear();
	std::cout <<"\n == CONNECTED COMPONENTS LEFT_SIDE_PATTERN ==\n" <<std::endl;
	ggl::LeftSidePattern ls(r);

	std::cout <<"--> connectedComponents( ls ) = "; std::cout.flush();
	std::cout	<<sgm::Graph_Interface::connectedComponents( ls, label ) 
				<<std::endl;
	
	std::cout <<"--> component labeling :" <<std::endl;
	for (size_t i=0; i<label.size(); i++)
		std::cout <<" " <<i <<" = " <<label[i] <<std::endl;


	std::cout <<"--> ls.getComponentLabeling() : " <<std::endl;
	for (size_t i=0; i<ls.getComponentLabeling().size(); i++)
		std::cout <<" " <<i <<" = " <<ls.getComponentLabeling()[i] <<std::endl;
	
	std::cout <<"--> ls.getFirstOfEachComponent() : "; std::cout.flush();
	const ggl::LeftSidePattern::IndexSet& idx = ls.getFirstOfEachComponent();
	for (	ggl::LeftSidePattern::IndexSet::const_iterator it = idx.begin();
			it != idx.end(); it++) {
		std::cout <<*it <<" "; std::cout.flush();
	}
	std::cout <<std::endl;
	
	label.clear();
	std::cout <<"\n == CONNECTED COMPONENTS RIGHT_SIDE_PATTERN ==\n" <<std::endl;
	ggl::RightSidePattern rs(r);

	std::cout <<"--> connectedComponents( rs ) = "; std::cout.flush();
	std::cout	<<sgm::Graph_Interface::connectedComponents( rs, label ) 
				<<std::endl;
	
	std::cout <<"--> component labeling :" <<std::endl;
	for (size_t i=0; i<label.size(); i++)
		std::cout <<" " <<i <<" = " <<label[i] <<std::endl;
		

}

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"                ggl::Rule  \n" 
				<<"==============================================\n" 
				<<std::endl;
	

	std::string ruleString;
	
	  // check rule 1
	{
		std::cout <<"\n == RULE ==\n" <<std::endl;
	
		ggl::Rule r(getRule_1(ruleString), "");
		
		std::cout <<"\n" <<ruleString <<"\n"  <<std::endl;
	
		printRule(r);
		printLeftSidePattern(r);
		printRightSidePattern(r);
		checkConnectedComponents(r);
	}
	
	  // check rule 2
	{
		std::cout <<"\n == RULE ==\n" <<std::endl;
	
		ggl::Rule r(getRule_2(ruleString), "rule2");
		
		std::cout <<"\n" <<ruleString <<"\n"  <<std::endl;
	
		printRule(r);
		printLeftSidePattern(r);
		printRightSidePattern(r);
		checkConnectedComponents(r);
	}
	
	  // check rule 3
	{
		std::cout <<"\n == RULE ==\n" <<std::endl;
	
		ggl::Rule r(getRule_3(ruleString), "rule3");
		
		std::cout <<"\n" <<ruleString <<"\n"  <<std::endl;
	
		printRule(r);
		printLeftSidePattern(r);
		printRightSidePattern(r);
		checkConnectedComponents(r);
	}
	
	  // check rule 4
	{
		std::cout <<"\n == RULE ==\n" <<std::endl;

		ggl::Rule r(getRule_4(ruleString),"rule4");

		std::cout <<"\n" <<ruleString <<"\n"  <<std::endl;

		printRule(r);
		printLeftSidePattern(r);
		printRightSidePattern(r);
		checkConnectedComponents(r);
	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
