
#include <iostream>
#include <sstream>

#include "ggl/Graph.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/Graph_GMLparser.hh"

#include "sgm/Graph_boost.hh"
#include "sgm/GraphScaffold.hh"

#include "utilPrintGraph_Interface.icc"

using namespace ggl;
using namespace sgm;

const std::string gml[] =  {
		"graph [ \n"
				" node [ id 0 label \"a\" ]\n"
				" node [ id 1 label \"b\" ]\n"
				" node [ id 2 label \"c\" ]\n"
				" node [ id 3 label \"d\" ]\n"
				" edge [ source 1 target 0 label \"-\" ]\n"
				" edge [ source 2 target 1 label \"-\" ]\n"
				" edge [ source 3 target 1 label \"-\" ]\n"
				"]\n"
		,
		"graph [ \n"
				" node [ id 0 label \"a\" ]\n"
				" node [ id 1 label \"b\" ]\n"
				" node [ id 2 label \"c\" ]\n"
				" node [ id 3 label \"d\" ]\n"
				" node [ id 4 label \"e\" ]\n"
				" node [ id 5 label \"f\" ]\n"
				" node [ id 6 label \"g\" ]\n"
				" node [ id 7 label \"h\" ]\n"
				" node [ id 8 label \"i\" ]\n"
				" node [ id 9 label \"j\" ]\n"
				" node [ id 10 label \"k\" ]\n"
				" node [ id 11 label \"l\" ]\n"
				" edge [ source 1 target 0 label \"-\" ]\n"
				" edge [ source 2 target 1 label \"-\" ]\n"
				" edge [ source 3 target 2 label \"-\" ]\n"
				" edge [ source 3 target 0 label \"-\" ]\n"
				" edge [ source 4 target 3 label \"-\" ]\n"
				" edge [ source 5 target 4 label \"-\" ]\n"
				" edge [ source 6 target 5 label \"-\" ]\n"
				" edge [ source 7 target 6 label \"-\" ]\n"
				" edge [ source 8 target 7 label \"-\" ]\n"
				" edge [ source 9 target 8 label \"-\" ]\n"
				" edge [ source 9 target 6 label \"-\" ]\n"
				" edge [ source 9 target 10 label \"-\" ]\n"
				" edge [ source 9 target 11 label \"-\" ]\n"
				" edge [ source 11 target 11 label \"-\" ]\n"
			"]\n"
		,
		""
};

 /*!
  * RingReporter that just counts the number of reported rings.
  */
class RR_Count : public sgm::RingReporter {
public:
	size_t ringCount;
	  //! construction
	RR_Count() : RingReporter(), ringCount(0) { reset(); }
	  //! destruction
	virtual ~RR_Count(){}

	void
	reset() {ringCount = 0;};

	void
	reportRing( const Graph_Interface& graph, const RingList & ringList )
	{
		ringCount++;
//		std::ostream& out = std::cout;
//		out <<"  * reported ring : ";
//		for (RingList::const_iterator it=ringList.begin(); it!=ringList.end(); it++) {
//			out <<*(it) <<" ";
//		}
//		out <<std::endl;
	}


};

int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"           sgm::GraphScaffold  \n"
				<<"==============================================\n" 
				<<std::endl;

	for (size_t i=0; !(gml[i].empty()); i++)
	{
		std::cout <<"\n ####### NEXT GRAPH ##################################\n"
				  <<std::endl;

		std::string graphString = gml[i];
		
		std::cout <<"\n graph to check:\n" <<graphString <<std::endl;
		
		std::cout <<"\n-->Graph_GMLparser::parseGraph( GML )" <<std::endl;
		
		std::pair<ggl::Graph, int> ret = ggl::Graph_GMLparser::parseGraph( graphString );
		
		if (ret.second >= 0 ) {
			std::cout <<"\n PARSING ERROR : at input position " <<ret.second
					<<" = '" <<graphString.at(ret.second) <<"' "
					<<" within \n'"
					<<graphString.substr(std::max(0,((int)ret.second)-10),20)
					<<"'\n"
					<<std::endl;
			continue;
		}

		sgm::Graph_boost<ggl::Graph> gi(ret.first);
		std::cout <<"\n result graph : " <<std::endl;
		print(gi);
		std::cout <<std::endl;

		  // create scaffold annotation object
		GraphScaffold scaffold;

		std::cout <<"\n-->GraphScaffold::getScaffoldAnnotation( graph )" <<std::endl;
		const GraphScaffold::ScaffoldAnnotation annotation
			= scaffold.getScaffoldAnnotation( gi );

		  // print annotation
		for (size_t i=0; i < annotation.size(); ++i) {
			switch (annotation.at(i)) {
				case GraphScaffold::GST_RING:
					std::cout <<"\t" <<i <<" = ring" <<std::endl;
					break;
				case GraphScaffold::GST_LINKER:
					std::cout <<"\t" <<i <<" = linker" <<std::endl;
					break;
				case GraphScaffold::GST_DANGLING:
					std::cout <<"\t" <<i <<" = dangling" <<std::endl;
					break;
				default:
					std::cout <<"\t" <<i <<" = unknown" <<std::endl;
					break;
			}
		}

	}


	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<std::endl;
	return 0;
}
