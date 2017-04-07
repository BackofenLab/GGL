
#include "sgm/HashMap.hh"

#if HAVE_UNORDERED_MAP > 0
	#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP > 0
	#include <tr1/unordered_map>
#elif HAVE_GNU_HASH_MAP > 0
	#include <ext/hash_map>
#else
	#include <map>
#endif


#include "ggl/Graph_GML_writer.hh"

namespace ggl
{

////////////////////////////////////////////////////////////////////////////////


void
Graph_GML_writer
::write( std::ostream& out
		, const Graph &graph
		, const bool withSpaces
		, const std::string* additionalGML )
{
	boost::property_map<Graph, PropNodeLabel>::const_type
		nodeLabel = boost::get( PropNodeLabel(), graph );
	boost::property_map<Graph, PropEdgeLabel>::const_type
		edgeLabel = boost::get( PropEdgeLabel(), graph );
	
	typedef
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map< Graph::vertex_descriptor, size_t>
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map< Graph::vertex_descriptor, size_t>
#else
		std::map< Graph::vertex_descriptor, size_t>
#endif
		Node2IndexMap;
		
	Node2IndexMap node2idx;
	
	  // open GML output
	out	<<"graph"
		<<(withSpaces?" ":"")
		<<"["
		<<(withSpaces?"\n":"");
	
	  // write nodes
	 Graph::vertex_iterator vi, v_end;
	boost::tie(vi,v_end) = boost::vertices(graph);
	size_t idx = 0;
	while (vi != v_end) {
		node2idx[*vi] = idx;
		
		out	<<(withSpaces?"  ":"")
			<<"node"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?" ":"")
			<<"id "
			<<idx
			<<" label \""
			<<nodeLabel[*vi]
			<<"\""
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
		idx++;
		vi++;
	}
	
	// write edges
	Graph::edge_iterator ei, e_end;
	boost::tie(ei,e_end) = boost::edges(graph);
	while (ei != e_end) {
		
		out	<<(withSpaces?"  ":"")
			<<"edge"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?" ":"")
			<<"source "
			<<node2idx[boost::source(*ei,graph)]
			<<" target "
			<<node2idx[boost::target(*ei,graph)]
			<<" label \""
			<<edgeLabel[*ei]
			<<"\""
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
		ei++;
	}
	
	 // add additional GML if present
	if (additionalGML != NULL)  {
		if (withSpaces) {
			  // just output additional GML and end with newline
			out <<*additionalGML <<"\n";
		} else {
			  // replace new lines with spaces
			std::string withoutNewline = *additionalGML;
			size_t newlinePos = withoutNewline.find('\n');
			while (newlinePos != std::string::npos) {
				withoutNewline[newlinePos] = ' ';
				newlinePos = withoutNewline.find('\n', newlinePos);
			}
			  // output stripped additional GML
			out <<withoutNewline;
		}
	}

	 // close GML output
	out <<"]"
		<<(withSpaces?"\n":"")
		;
	
}

////////////////////////////////////////////////////////////////////////////////


} // namespace ggl
