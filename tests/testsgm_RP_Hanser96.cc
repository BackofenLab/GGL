
#include <iostream>
#include <sstream>

#include "ggl/Graph.hh"
#include "ggl/Graph_GML_writer.hh"
#include "ggl/Graph_GMLparser.hh"

#include "sgm/Graph_boost.hh"
#include "sgm/RP_Hanser96.hh"

#include "dataTargetGraph_1.icc"
#include "dataTargetGraph_2.icc"
#include "dataTargetGraph_3.icc"

#include "utilPrintGraph_Interface.icc"

using namespace ggl;
using namespace sgm;

const std::string gml[] =  {
		"graph [ \n"
				" node [ id 0 label \"a\" ]\n"
				" node [ id 1 label \"b\" ]\n"
				" node [ id 2 label \"c\" ]\n"
				" node [ id 3 label \"d\" ]\n"
				" node [ id 4 label \"e\" ]\n"
				" node [ id 5 label \"f\" ]\n"
				" node [ id 6 label \"g\" ]\n"
				" node [ id 7 label \"h\" ]\n"
				" edge [ source 0 target 1 label \"-\" ]\n"
				" edge [ source 0 target 2 label \"-\" ]\n"
				" edge [ source 0 target 7 label \"-\" ]\n"
				" edge [ source 1 target 2 label \"-\" ]\n"
				" edge [ source 1 target 3 label \"-\" ]\n"
				" edge [ source 3 target 4 label \"-\" ]\n"
				" edge [ source 3 target 5 label \"-\" ]\n"
				" edge [ source 3 target 6 label \"-\" ]\n"
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
				" edge [ source 0 target 3 label \"-\" ]\n"
				" edge [ source 0 target 6 label \"-\" ]\n"
				" edge [ source 0 target 4 label \"-\" ]\n"
				" edge [ source 1 target 4 label \"-\" ]\n"
				" edge [ source 1 target 5 label \"-\" ]\n"
				" edge [ source 1 target 7 label \"-\" ]\n"
				" edge [ source 2 target 3 label \"-\" ]\n"
				" edge [ source 2 target 5 label \"-\" ]\n"
				" edge [ source 2 target 8 label \"-\" ]\n"
				" edge [ source 3 target 6 label \"-\" ]\n"
				" edge [ source 3 target 8 label \"-\" ]\n"
				" edge [ source 4 target 6 label \"-\" ]\n"
				" edge [ source 4 target 7 label \"-\" ]\n"
				" edge [ source 5 target 7 label \"-\" ]\n"
				" edge [ source 5 target 8 label \"-\" ]\n"
				" edge [ source 6 target 7 label \"-\" ]\n"
				" edge [ source 6 target 8 label \"-\" ]\n"
				" edge [ source 7 target 8 label \"-\" ]\n"
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
				" node [ id 12 label \"m\" ]\n"
				" node [ id 13 label \"n\" ]\n"
				" node [ id 14 label \"o\" ]\n"
				" node [ id 15 label \"p\" ]\n"
				" node [ id 16 label \"q\" ]\n"
				" edge [ source 0 target 1 label \"-\" ]\n"
				" edge [ source 0 target 5 label \"-\" ]\n"
				" edge [ source 1 target 2 label \"-\" ]\n"
				" edge [ source 2 target 3 label \"-\" ]\n"
				" edge [ source 2 target 6 label \"-\" ]\n"
				" edge [ source 3 target 4 label \"-\" ]\n"
				" edge [ source 3 target 9 label \"-\" ]\n"
				" edge [ source 4 target 5 label \"-\" ]\n"
				" edge [ source 6 target 7 label \"-\" ]\n"
				" edge [ source 6 target 10 label \"-\" ]\n"
				" edge [ source 7 target 8 label \"-\" ]\n"
				" edge [ source 7 target 13 label \"-\" ]\n"
				" edge [ source 8 target 9 label \"-\" ]\n"
				" edge [ source 10 target 11 label \"-\" ]\n"
				" edge [ source 11 target 12 label \"-\" ]\n"
				" edge [ source 12 target 13 label \"-\" ]\n"
				" edge [ source 12 target 14 label \"-\" ]\n"
				" edge [ source 13 target 16 label \"-\" ]\n"
				" edge [ source 14 target 15 label \"-\" ]\n"
				" edge [ source 15 target 16 label \"-\" ]\n"
			"]\n"
		,
		"graph [\n"
				"  node [ id 1 label \"1\" ] \n"
				"  node [ id 2 label \"2\" ] \n"
				"  node [ id 3 label \"3\" ] \n"
				"  node [ id 4 label \"4\" ] \n"
				"  node [ id 5 label \"5\" ] \n"
				"  node [ id 6 label \"6\" ] \n"
				"  node [ id 7 label \"7\" ] \n"
				"  node [ id 8 label \"8\" ] \n"
				"  edge [ source 1 target 2 label \"-\" ] \n"
				"  edge [ source 1 target 3 label \"-\" ] \n"
				"  edge [ source 1 target 5 label \"-\" ] \n"
				"  edge [ source 2 target 4 label \"-\" ] \n"
				"  edge [ source 2 target 6 label \"-\" ] \n"
				"  edge [ source 3 target 4 label \"-\" ] \n"
				"  edge [ source 3 target 7 label \"-\" ] \n"
				"  edge [ source 4 target 8 label \"-\" ] \n"
				"  edge [ source 5 target 6 label \"-\" ] \n"
				"  edge [ source 5 target 7 label \"-\" ] \n"
				"  edge [ source 6 target 8 label \"-\" ] \n"
				"  edge [ source 7 target 8 label \"-\" ]\n"
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
				<<"           sgm::RP_Hanser96  \n"
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

		  // create ring perception object
		RP_Hanser96 hanser96;
		for (size_t maxRingSize = gi.getNodeNumber(); maxRingSize>0; --maxRingSize) {
			  // set up reporter
			RR_Count reporter;
			reporter.reset();

			std::cout <<"\n-->RP_Hanser96::findRings( graph, RR_Count, "<<maxRingSize<<" )" <<std::endl;
			const size_t ringCount = hanser96.findRings( gi, reporter, maxRingSize );

			std::cout	<<" + returned ring number = "<<ringCount <<std::endl
						<<" + reported ring number = "<<reporter.ringCount <<std::endl;
		}

		RP_Hanser96::BondSet ringBonds;
		const size_t ringBondCount = hanser96.findRingBonds( gi, ringBonds);
		std::cout	<<" + ring bond number     = "<<ringBondCount  <<std::endl;

	}


	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<std::endl;
	return 0;
}
