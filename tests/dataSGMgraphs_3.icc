#ifndef DATASGMGRAPHS_3_ICC_
#define DATASGMGRAPHS_3_ICC_

#ifndef DATASGMGRAPHS_GRAPH_CLASS
#define DATASGMGRAPHS_GRAPH_CLASS

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
	
	//! The properties available for the nodes of a graph
	typedef	boost::property<	boost::vertex_index_t, size_t
			, boost::property<	boost::vertex_name_t,  std::string
				> >
					MyGraph_NodeProperties;
	
	//! The properties available for the edges of a graph
	typedef	boost::property<	boost::edge_name_t, std::string 
				> 
					MyGraph_EdgeProperties;
	
	//! The definition of undirected graphs 
	typedef boost::adjacency_list<
					boost::vecS,      				// store edges
					boost::vecS,       				// store vertices
					boost::undirectedS,				// is an undirected graph
					MyGraph_NodeProperties,  		// (atom symbols etc) 
					MyGraph_EdgeProperties   		// (edge symbols etc)
				> 
					MyGraph;


#endif // DATASGMGRAPHS_GRAPH_CLASS

	MyGraph
	getPattern_3 (std::string& graphString) {
		
		  // the target graph to fill
		MyGraph t;
		
		MyGraph::vertex_descriptor v;

		boost::property_map< MyGraph, boost::vertex_name_t >::type 
			nodeLabel = boost::get( boost::vertex_name_t(), t );
		
		MyGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type 
			edgeLabel = boost::get( boost::edge_name_t(), t );
		
		graphString = 
"\n=============================================\n\n"
"  0(C)- \n"
"    |  |\n"
"    1 - \n"
"\n=============================================\n";

		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
			
		// edge 0 - 0
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(0,t), t).first;
		edgeLabel[e] = "-1-";
		
		return t;
	}


	MyGraph
	getPattern_3_WC (std::string& graphString) {

		  // the target graph to fill
		MyGraph t;

		MyGraph::vertex_descriptor v;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), t );

		MyGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type
			edgeLabel = boost::get( boost::edge_name_t(), t );

		graphString =
"\n=============================================\n\n"
"  0(*)- \n"
"    |  |\n"
"    * - \n"
"\n=============================================\n";

		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "*";

		// edge 0 - 0
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(0,t), t).first;
		edgeLabel[e] = "*";

		return t;
	}

	MyGraph
	getTarget_3 (std::string& graphString) {
		
		  // the target graph to fill
		MyGraph t;
		
		MyGraph::vertex_descriptor v;

		boost::property_map< MyGraph, boost::vertex_name_t >::type 
			nodeLabel = boost::get( boost::vertex_name_t(), t );
		
		MyGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type 
			edgeLabel = boost::get( boost::edge_name_t(), t );
		
		graphString = 
				"\n=============================================\n\n"
				"  0(C)-1\n"
				"    | \\/\n"
				"    1 \n"
				"    | \n"
				"  1(N)-1\n"
				"    | \\/\n"
				"    2 \n"
				"    | \n"
				"  2(N)-2\n"
				"      \\/\n"
				"\n=============================================\n";

		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 1
		v = boost::add_vertex(t);
		nodeLabel[v] = "N";
		// node 2
		v = boost::add_vertex(t);
		nodeLabel[v] = "N";
			
		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 1 - 2
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(2,t), t).first;
		edgeLabel[e] = "-2-";
		// loop 0
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(0,t), t).first;
		edgeLabel[e] = "-1-";
		// loop 1
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "-1-";
		// loop 2
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(2,t), t).first;
		edgeLabel[e] = "-2-";
		
		return t;
	}



#endif /*DATASGMGRAPHS_3_ICC_*/
