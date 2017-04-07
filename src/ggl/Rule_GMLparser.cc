
#include "ggl/Rule_GMLparser.hh"
#include "ggl/Rule_GML_grammar.hh"
#include "ggl/Graph.hh"
#include "ggl/Graph_GMLparser.hh"


namespace ggl {

//##############################################################################


	std::pair< Rule, int >
	Rule_GMLparser
    ::parseRule( const std::string & GML_string ) throw (Rule_GML_error)
    {
		  // forward call
		return Rule_GML_grammar::parseRule( GML_string );
    }

//##############################################################################


	std::pair< Rule, int >
	Rule_GMLparser
    ::parseCompactRule( const std::string & GML_string ) throw (Rule_GML_error)
    {
		  // forward call for parsing
		std::pair< Graph, int > graphParse;
		try {
			graphParse = Graph_GMLparser::parseGraph( GML_string );
		} catch ( std::exception ex ) {
			throw new Rule_GML_error(ex.what());
		}
		ggl::Graph & graph = graphParse.first;
		// Access to the node label property_map of graph to set node labels
		boost::property_map<Graph, PropNodeLabel>::type graphNodeLabel(boost::get( PropNodeLabel(), graph ));
		// Access to the node index property_map of graph to set node labels
		boost::property_map<Graph, PropNodeIndex>::const_type graphNodeIndex(boost::get( PropNodeIndex(), graph ));
		// Access to the edge label property_map of graph to set edge labels
		boost::property_map<Graph, PropEdgeLabel>::type graphEdgeLabel(boost::get( PropEdgeLabel(), graph ));


		Rule::CoreGraph core;

		  // parsing error handling
		if (graphParse.second != -1) {
			return std::pair< Rule, int>(Rule(core,"",NULL),graphParse.second);
		}

		// Access to the node label property_map of core to set node labels
		boost::property_map<Rule::CoreGraph, Rule::NodeLabelProperty>::type coreNodeLabel(boost::get( Rule::NodeLabelProperty(), core ));
		// Access to the right node label property_map of core to set changing node labels
		boost::property_map<Rule::CoreGraph, Rule::NodeRightLabelProperty>::type coreNodeLabelRight(boost::get( Rule::NodeRightLabelProperty(), core ));
		// Access to the node Rule context property_map of core
		boost::property_map<Rule::CoreGraph, Rule::NodeContextProperty>::type coreNodeContext(boost::get( Rule::NodeContextProperty(), core ));
		// Access to the edge label property_map of core to set edge labels
		boost::property_map<Rule::CoreGraph, Rule::EdgeLabelProperty>::type coreEdgeLabel(boost::get( Rule::EdgeLabelProperty(), core ));
		// Access to the edge Rule context property_map of core
		boost::property_map<Rule::CoreGraph, Rule::EdgeContextProperty>::type coreEdgeContext(boost::get( Rule::EdgeContextProperty(), core ));


		// copy nodes
		Rule::CoreGraph::vertex_descriptor coreNode;
		Graph::vertex_iterator vi, v_end;
		for (boost::tie(vi,v_end) = boost::vertices(graph); vi != v_end;++vi) {

			  // create new node within core graph
			coreNode = boost::add_vertex( core );

			  // set node label depending on separator
			if ( graphNodeLabel[*vi].find('|') != std::string::npos ) {
				size_t splitPos = graphNodeLabel[*vi].find('|');
				std::string firstLabel = graphNodeLabel[*vi].substr(0,splitPos);
				std::string newLabel = graphNodeLabel[*vi].substr(splitPos+1);
				if (firstLabel.empty()) {
					  // copy label
					coreNodeLabel[coreNode] = newLabel;
					  // set to right side only
					coreNodeContext[coreNode] = Rule::RULE_RIGHT_SIDE;
				} else if (newLabel.empty()) {
					  // copy label
					coreNodeLabel[coreNode] = firstLabel;
					  // set to left side only
					coreNodeContext[coreNode] = Rule::RULE_LEFT_SIDE;
				} else {
					  // copy label
					coreNodeLabel[coreNode] = firstLabel;
					coreNodeLabelRight[coreNode] = newLabel;
					  // set to label change
					coreNodeContext[coreNode] = Rule::RULE_LABEL_CHANGE;
				}
			} else {
				  // copy label
				coreNodeLabel[coreNode] = graphNodeLabel[*vi];
				  // set to context
				coreNodeContext[coreNode] = Rule::RULE_CONTEXT;
			}
		}

		// copy edges
		Rule::CoreGraph::edge_descriptor coreEdge;
		Graph::edge_iterator ei, e_end;
		for (boost::tie(ei,e_end) = boost::edges(graph); ei != e_end; ++ei) {
			  // add new edge
			coreEdge = boost::add_edge(
					boost::vertex( graphNodeIndex[boost::source(*ei,graph)], core),
					boost::vertex( graphNodeIndex[boost::target(*ei,graph)], core),
					core
					).first;

			  // set edge label depending on separator
			if ( graphEdgeLabel[*ei].find('|') != std::string::npos ) {
				size_t splitPos = graphEdgeLabel[*ei].find('|');
				std::string firstLabel = graphEdgeLabel[*ei].substr(0,splitPos);
				std::string newLabel = graphEdgeLabel[*ei].substr(splitPos+1);
				if (firstLabel.empty()) {
					  // copy label
					coreEdgeLabel[coreEdge] = newLabel;
					  // set to right side only
					coreEdgeContext[coreEdge] = Rule::RULE_RIGHT_SIDE;
				} else if (newLabel.empty()) {
					  // copy label
					coreEdgeLabel[coreEdge] = firstLabel;
					  // set to left side only
					coreEdgeContext[coreEdge] = Rule::RULE_LEFT_SIDE;
				} else {
					  // copy label
					coreEdgeLabel[coreEdge] = firstLabel;
					  // set to label change
					coreEdgeContext[coreEdge] = Rule::RULE_LEFT_SIDE;
					  // add new edge for right side label
					coreEdge = boost::add_edge(
							boost::vertex( graphNodeIndex[boost::source(*ei,graph)], core),
							boost::vertex( graphNodeIndex[boost::target(*ei,graph)], core),
							core
							).first;
					  // copy label
					coreEdgeLabel[coreEdge] = newLabel;
					  // set to label change
					coreEdgeContext[coreEdge] = Rule::RULE_RIGHT_SIDE;
				}
			} else {
				  // copy label
				coreEdgeLabel[coreEdge] = graphEdgeLabel[*ei];
				  // set to context
				coreEdgeContext[coreEdge] = Rule::RULE_CONTEXT;
			}

		}

		return std::pair< Rule, int>(Rule(core,"",NULL),-1);
    }

//##############################################################################

} // namespace ggl

