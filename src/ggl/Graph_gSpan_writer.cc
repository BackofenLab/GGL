
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


#include "ggl/Graph_gSpan_writer.hh"

namespace ggl
{

////////////////////////////////////////////////////////////////////////////////


void
Graph_gSpan_writer
::write( std::ostream& out, const Graph &graph)
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
	
	  // open gSpan output
	out	<<"t\n";
	
	  // write nodes
	 Graph::vertex_iterator vi, v_end;
	boost::tie(vi,v_end) = boost::vertices(graph);
	size_t idx = 0;
	while (vi != v_end) {
		node2idx[*vi] = idx;
		  // create vertex entry
		out	<<"v " <<idx <<" " <<nodeLabel[*vi] <<"\n";
		  // increase counter
		idx++;
		vi++;
	}
	
	// write edges
	Graph::edge_iterator ei, e_end;
	boost::tie(ei,e_end) = boost::edges(graph);
	while (ei != e_end) {
		
		out	<<"e "
			<<node2idx[boost::source(*ei,graph)]
			<<" "
			<<node2idx[boost::target(*ei,graph)]
			<<" "
			<<edgeLabel[*ei]
			<<"\n"
			;
		ei++;
	}
	
}

////////////////////////////////////////////////////////////////////////////////


} // namespace ggl
