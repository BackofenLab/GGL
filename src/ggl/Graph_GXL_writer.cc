
#include "ggl/Graph_GXL_writer.hh"

namespace ggl
{

////////////////////////////////////////////////////////////////////////////////


void
Graph_GXL_writer
::write( std::ostream& out, const Graph &graph )
{
	boost::property_map<Graph, PropNodeLabel>::const_type
		nodeLabel = boost::get( PropNodeLabel(), graph );
	boost::property_map<Graph, PropEdgeLabel>::const_type
		edgeLabel = boost::get( PropEdgeLabel(), graph );
	boost::property_map<Graph, PropNodeIndex>::const_type
		nodeIndex = boost::get( PropNodeIndex(), graph );
	
		
	  // open GXL output
	out	<<	"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n"
			"<gxl>\n"
			" <graph id=\"GGL-graph\" edgeids=\"false\" edgemode=\"undirected\" hypergraph=\"false\">\n"
			;
	
	  // write nodes
	Graph::vertex_iterator vi, v_end;

	for (boost::tie(vi,v_end) = boost::vertices(graph); vi != v_end; ++vi) {
		
		out	<<	"  <node id=\""
			<<	nodeIndex[*vi]
			<<	"\">"
				" <attr name=\"label\">"
				" <string>"
			<<	nodeLabel[*vi]
			<<	"</string>"
				" </attr>"
				" </node>\n"
			;
	}
	
	// write edges
	Graph::edge_iterator ei, e_end;

	for (boost::tie(ei,e_end) = boost::edges(graph); ei != e_end; ++ei) {
		
		out	<<	"  <edge to=\""
			<<	nodeIndex[boost::target(*ei,graph)]
			<<	"\" from=\""
			<<	nodeIndex[boost::source(*ei,graph)]
			<<	"\">"
				" <attr name=\"label\">"
				" <string>"
			<<	edgeLabel[*ei]
			<<	"</string>"
				" </attr>"
				" </edge>\n"
			;
	}
	
	 // close GXL output
	out <<	" </graph>\n"
			"</gxl>\n"
			;
	
}

////////////////////////////////////////////////////////////////////////////////


} // namespace ggl
