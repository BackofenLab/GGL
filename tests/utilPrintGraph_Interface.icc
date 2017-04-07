
#include <iostream>
#include <iomanip>

#include "sgm/Graph_Interface.hh"

	void
	print( const sgm::Graph_Interface& g ) 
	{
		for (size_t i=0; i<g.getNodeNumber(); i++) {
			std::cout <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i);
					e != g.getOutEdgesEnd(i); ++e )
			{
				std::cout <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			std::cout <<" |\n";
		}

	}


