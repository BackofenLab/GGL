#ifndef DATARULE_1_ICC_
#define DATARULE_1_ICC_

#include <ggl/Rule.hh>


	ggl::Rule::CoreGraph
	getRule_1 (std::string & ruleString) {
		
		  // fill rule graph
		ggl::Rule::CoreGraph g;
		{
		ggl::Rule::CoreGraph::vertex_descriptor v;
		  // get level property class
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::NodeContextProperty>::type
			nodeLevel = boost::get( ggl::Rule::NodeContextProperty(), g );
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::NodeLabelProperty >::type
			nodeLabel = boost::get( ggl::Rule::NodeLabelProperty(), g );
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::NodeRightLabelProperty >::type
			nodeRightLabel = boost::get( ggl::Rule::NodeRightLabelProperty(), g );
		
		ggl::Rule::CoreGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::EdgeContextProperty>::type
			edgeLevel = boost::get( ggl::Rule::EdgeContextProperty(), g );
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::EdgeLabelProperty >::type
			edgeLabel = boost::get( ggl::Rule::EdgeLabelProperty(), g );
	
		ruleString = 
"\n=============================================\n\n\
  0(A) -1- 1(B)               0(E) -1- 1(B)\n\
   |                           |        |\n\
   2               ==>         4        3\n\
   |                           |        |\n\
  2(C)                        2(C)     3(D)\n\
\n=============================================\n";
	
		// node 0
		v = boost::add_vertex(g);
		nodeLabel[v] = "A";
		nodeRightLabel[v] = "E";
		nodeLevel[v] = ggl::Rule::RULE_LABEL_CHANGE;
		// node 1
		v = boost::add_vertex(g);
		nodeLabel[v] = "B";
		nodeLevel[v] = ggl::Rule::RULE_CONTEXT;
		// node 2
		v = boost::add_vertex(g);
		nodeLabel[v] = "C";
		nodeLevel[v] = ggl::Rule::RULE_CONTEXT;
		// node 3
		v = boost::add_vertex(g);
		nodeLabel[v] = "D";
		nodeLevel[v] = ggl::Rule::RULE_RIGHT_SIDE;
		
		//	0(A) -1- 1(B)				0(E) -1- 1(B)
		//	 |							 |        |
		//	 2					==>		 4        3
		//	 |							 |        |
		//	2(C)						2(C)     3(D)
	
		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(1,g), g).first;
		edgeLabel[e] = "-1-";
		edgeLevel[e] = ggl::Rule::RULE_CONTEXT;
		// edge 0 - 2
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(2,g), g).first;
		edgeLabel[e] = "-2-";
		edgeLevel[e] = ggl::Rule::RULE_LEFT_SIDE;
		// edge 0 - 2
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(2,g), g).first;
		edgeLabel[e] = "-4-";
		edgeLevel[e] = ggl::Rule::RULE_RIGHT_SIDE;
		// edge 1 - 3
		e = boost::add_edge(boost::vertex(1,g), boost::vertex(3,g), g).first;
		edgeLabel[e] = "-3-";
		edgeLevel[e] = ggl::Rule::RULE_RIGHT_SIDE;
		
		}
		
		return g;
	}

	ggl::Rule::CoreGraph
	getRule_1_WC (std::string & ruleString) {

		  // fill rule graph
		ggl::Rule::CoreGraph g;
		{
		ggl::Rule::CoreGraph::vertex_descriptor v;
		  // get level property class
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::NodeContextProperty>::type
			nodeLevel = boost::get( ggl::Rule::NodeContextProperty(), g );
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::NodeLabelProperty >::type
			nodeLabel = boost::get( ggl::Rule::NodeLabelProperty(), g );
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::NodeRightLabelProperty >::type
			nodeRightLabel = boost::get( ggl::Rule::NodeRightLabelProperty(), g );

		ggl::Rule::CoreGraph::edge_descriptor e;
		  // get level property class
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::EdgeContextProperty>::type
			edgeLevel = boost::get( ggl::Rule::EdgeContextProperty(), g );
		boost::property_map< ggl::Rule::CoreGraph, ggl::Rule::EdgeLabelProperty >::type
			edgeLabel = boost::get( ggl::Rule::EdgeLabelProperty(), g );

		ruleString =
"\n=============================================\n\n\
  0(*) -1- 1(B)               0(*) -1- 1(E)\n\
   |                           |        |\n\
   2               ==>         4        3\n\
   |                           |        |\n\
  2(C)                        2(C)     3(D)\n\
\n=============================================\n";

		// node 0
		v = boost::add_vertex(g);
		nodeLabel[v] = "*";
		nodeLevel[v] = ggl::Rule::RULE_CONTEXT;
		// node 1
		v = boost::add_vertex(g);
		nodeLabel[v] = "B";
		nodeRightLabel[v] = "E";
		nodeLevel[v] = ggl::Rule::RULE_LABEL_CHANGE;
		// node 2
		v = boost::add_vertex(g);
		nodeLabel[v] = "C";
		nodeLevel[v] = ggl::Rule::RULE_CONTEXT;
		// node 3
		v = boost::add_vertex(g);
		nodeLabel[v] = "D";
		nodeLevel[v] = ggl::Rule::RULE_RIGHT_SIDE;

		//	0(A) -1- 1(B)				0(E) -1- 1(B)
		//	 |							 |        |
		//	 2					==>		 4        3
		//	 |							 |        |
		//	2(C)						2(C)     3(D)

		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(1,g), g).first;
		edgeLabel[e] = "-1-";
		edgeLevel[e] = ggl::Rule::RULE_CONTEXT;
		// edge 0 - 2
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(2,g), g).first;
		edgeLabel[e] = "-2-";
		edgeLevel[e] = ggl::Rule::RULE_LEFT_SIDE;
		// edge 0 - 2
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(2,g), g).first;
		edgeLabel[e] = "-4-";
		edgeLevel[e] = ggl::Rule::RULE_RIGHT_SIDE;
		// edge 1 - 3
		e = boost::add_edge(boost::vertex(1,g), boost::vertex(3,g), g).first;
		edgeLabel[e] = "-3-";
		edgeLevel[e] = ggl::Rule::RULE_RIGHT_SIDE;

		}

		return g;
	}

#endif /*DATARULE_1_ICC_*/
