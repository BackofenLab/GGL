#ifndef DATARULE_3_ICC_
#define DATARULE_3_ICC_

#include <ggl/Rule.hh>


	ggl::Rule::CoreGraph
	getRule_3 (std::string & ruleString) {
		
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
  2(D) -1- 1(B)               2(C)     1(B) \n\
                                        |   \n\
                    ==>                 3   \n\
                                        |   \n\
           0(A)                        0(A) \n\
\n=============================================\n";
	
		// node 0
		v = boost::add_vertex(g);
		nodeLabel[v] = "A";
		nodeLevel[v] = ggl::Rule::RULE_CONTEXT;
		// node 1
		v = boost::add_vertex(g);
		nodeLabel[v] = "B";
		nodeLevel[v] = ggl::Rule::RULE_CONTEXT;
		// node 2
		v = boost::add_vertex(g);
		nodeLabel[v] = "D";
		nodeRightLabel[v] = "C";
		nodeLevel[v] = ggl::Rule::RULE_LABEL_CHANGE;
		
	
		// edge 1 - 2
		e = boost::add_edge(boost::vertex(1,g), boost::vertex(2,g), g).first;
		edgeLabel[e] = "-1-";
		edgeLevel[e] = ggl::Rule::RULE_LEFT_SIDE;
		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,g), boost::vertex(1,g), g).first;
		edgeLabel[e] = "-3-";
		edgeLevel[e] = ggl::Rule::RULE_RIGHT_SIDE;
		
		}
		
		return g;
	}

#endif /*DATARULE_3_ICC_*/
