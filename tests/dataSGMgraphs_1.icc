#ifndef DATASGMGRAPHS_1_ICC_
#define DATASGMGRAPHS_1_ICC_

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
	getPattern_1 (std::string& graphString) {
		
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
"\n=============================================\n\n\
  0(C)     1(C) \n\
   |      /     \n\
   1    1       \n\
   |  /         \n\
  3(N) -1- 2(C) \n\
\n=============================================\n";

		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 1
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 2
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 3
		v = boost::add_vertex(t);
		nodeLabel[v] = "N";
			
		// edge 0 - 3
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 1 - 3
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";
		
		return t;
	}


	MyGraph
	getPattern_1_WC (std::string& graphString) {

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
"\n=============================================\n\n\
  0(C)     1(*) \n\
   |      /     \n\
   1    1       \n\
   |  /         \n\
  3(N) -1- 2(C) \n\
\n=============================================\n";

		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 1
		v = boost::add_vertex(t);
		nodeLabel[v] = "*";
		// node 2
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 3
		v = boost::add_vertex(t);
		nodeLabel[v] = "N";

		// edge 0 - 3
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 1 - 3
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";

		return t;
	}

	MyGraph
	getTarget_1 (std::string& graphString) {
		
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
"\n=============================================\n\n\
  5(C) -1- 4(C) -2- 3(C)          \n\
   |                 |            \n\
   1                 1            \n\
   |                 |            \n\
  0(N) -1- 1(C) -1- 2(N)          \n\
   |                 |            \n\
   1                 1            \n\
   |                 |            \n\
  6(C)     8(O) -2- 7(C) -1- 9(C) \n\
\n=============================================\n";

		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "N";
		// node 1
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 2
		v = boost::add_vertex(t);
		nodeLabel[v] = "N";
		// node 3
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 4
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 5
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 6
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 7
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 8
		v = boost::add_vertex(t);
		nodeLabel[v] = "O";
		// node 9
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
			
		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 0 - 5
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(5,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 0 - 6
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(6,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 1 - 2
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(2,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 2 - 7
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(7,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 3 - 4
		e = boost::add_edge(boost::vertex(3,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 4 - 5
		e = boost::add_edge(boost::vertex(4,t), boost::vertex(5,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 7 - 8
		e = boost::add_edge(boost::vertex(7,t), boost::vertex(8,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 7 - 9
		e = boost::add_edge(boost::vertex(7,t), boost::vertex(9,t), t).first;
		edgeLabel[e] = "-1-";
		
		return t;
	}



#endif /*DATASGMGRAPHS_1_ICC_*/
