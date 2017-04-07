

#include <iomanip>
#include "ggl/Graph.hh"
#include "ggl/GS_stream.hh"

namespace ggl {

////////////////////////////////////////////////////////////////////////////////


	void
	GS_stream
	::add( const Graph & graph )
	{
		  // create graph representation that is easier to handle
		Graph_boost g(graph);
		  // write graph to stream
		for (size_t i=0; i<g.getNodeNumber(); ++i ) {
			out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
					eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
			{
				out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			out <<" |\n";
		}
		out <<std::endl;
	}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

} // namespace ggl

