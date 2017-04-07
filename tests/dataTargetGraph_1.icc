#ifndef DATATARGETGRAPH_1_ICC_
#define DATATARGETGRAPH_1_ICC_

#include <ggl/Graph.hh>

	ggl::Graph
	getTarget_1 (std::string& graphString) {
		
		  // the target graph to fill
		ggl::Graph t;
		
		ggl::Graph::vertex_descriptor v;
		boost::property_map< ggl::Graph, ggl::PropNodeLabel >::type
			nodeLabel = boost::get( ggl::PropNodeLabel(), t );
		
		ggl::Graph::edge_descriptor e;
		boost::property_map< ggl::Graph, ggl::PropEdgeLabel >::type
			edgeLabel = boost::get( ggl::PropEdgeLabel(), t );
		
		graphString = 
"\n=============================================\n\n\
  0(A) -1- 1(B) -1- 2(A)\n\
   |        |      /    \n\
   2        2    2      \n\
   |        |  /        \n\
  3(C) -1- 4(C)         \n\
\n=============================================\n";
	
		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "A";
		// node 1
		v = boost::add_vertex(t);
		nodeLabel[v] = "B";
		// node 2
		v = boost::add_vertex(t);
		nodeLabel[v] = "A";
		// node 3
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 4
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
	
		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 1 - 2
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(2,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 0 - 3
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 1 - 4
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 2 - 4
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 3 - 4
		e = boost::add_edge(boost::vertex(3,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-1-";
		
		return t;
	}



#endif /*DATATARGETGRAPH_1_ICC_*/
