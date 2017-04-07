#ifndef DATATARGETGRAPH_3_ICC_
#define DATATARGETGRAPH_3_ICC_

#include <ggl/Graph.hh>

	ggl::Graph
	getTarget_3 (std::string& graphString) {
		
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
  0(B) -1- 1(D) -1- 2(B) \n\
   | \\      |      / |  \n\
   2   3    2    3   1  \n\
   |     \\  |  /     |  \n\
  3(C) -1- 4(C) -0- 5(A) \n\
\n=============================================\n";
	
		// node 0
		v = boost::add_vertex(t);
		nodeLabel[v] = "B";
		// node 1
		v = boost::add_vertex(t);
		nodeLabel[v] = "D";
		// node 2
		v = boost::add_vertex(t);
		nodeLabel[v] = "B";
		// node 3
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 4
		v = boost::add_vertex(t);
		nodeLabel[v] = "C";
		// node 5
		v = boost::add_vertex(t);
		nodeLabel[v] = "A";

//		"\n=====================
//		  0(B) -1- 1(D) -1- 2(B)
//		   | \      |      /|   
//		   2   3    2    3  1   
//		   |     \  |  /    |   
//		  3(C) -1- 4(C) -0- 5(A)
//		\n======================
			
		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(1,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 0 - 3
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(3,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 0 - 4
		e = boost::add_edge(boost::vertex(0,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-3-";
		// edge 1 - 2
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(2,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 1 - 4
		e = boost::add_edge(boost::vertex(1,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-2-";
		// edge 2 - 4
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-3-";
		// edge 2 - 5
		e = boost::add_edge(boost::vertex(2,t), boost::vertex(5,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 3 - 4
		e = boost::add_edge(boost::vertex(3,t), boost::vertex(4,t), t).first;
		edgeLabel[e] = "-1-";
		// edge 4 - 5
		e = boost::add_edge(boost::vertex(4,t), boost::vertex(5,t), t).first;
		edgeLabel[e] = "-0-";
		
		return t;
	}



#endif /*DATATARGETGRAPH_3_ICC_*/
