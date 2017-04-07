
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

#include <sstream>
#include <stdexcept>

#include "ggl/Rule_GML_writer.hh"

namespace ggl
{

////////////////////////////////////////////////////////////////////////////////


void
Rule_GML_writer
::write( std::ostream& out, const Rule &rule, const bool withSpaces)
{
	  // access to core graph encoding the rule
	const Rule::CoreGraph & graph = rule.getCore();

	// TODO add constraints GML printing
	// TODO add copy and paste GML
	
	  // open GML output
	out	<<"rule"
		<<(withSpaces?" ":"")
		<<"["
		<<(withSpaces?"\n":"");
	
	  // write rule ID
	out	<<(withSpaces?" ":"")
		<<"ruleID \""
		<<rule.getID()
		<<"\""
		<<(withSpaces?"\n":"");

	  // print context
	std::string context = getContextGML( graph, Rule::RULE_CONTEXT, withSpaces );
	if (!context.empty()) {
		out <<(withSpaces?" ":"")
			<<"context"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?"\n":"")
			<<context
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
	}

	  // print left side
	context = getContextGML( graph, Rule::RULE_LEFT_SIDE, withSpaces );
	if (!context.empty()) {
		out <<(withSpaces?" ":"")
			<<"left"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?"\n":"")
			<<context
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
	}

	  // print right side
	context = getContextGML( graph, Rule::RULE_RIGHT_SIDE, withSpaces );
	if (!context.empty()) {
		out <<(withSpaces?" ":"")
			<<"right"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?"\n":"")
			<<context
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
	}

	 // close GML output
	out <<"]"
		<<(withSpaces?"\n":"")
		;

}

////////////////////////////////////////////////////////////////////////////////


void
Rule_GML_writer
::writeCompact( std::ostream& out, const Rule &rule, const bool withSpaces)
{
	  // access to core graph encoding the rule
	const Rule::CoreGraph & graph = rule.getCore();

	// TODO add constraints GML printing
	// TODO add copy and paste GML

	  // open GML output
	out	<<"graph"
		<<(withSpaces?" ":"")
		<<"["
		<<(withSpaces?"\n":"");

	boost::property_map<Rule::CoreGraph, Rule::NodeContextProperty>::const_type
		nodeContext = boost::get( Rule::NodeContextProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::NodeLabelProperty>::const_type
		nodeLabel = boost::get( Rule::NodeLabelProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::NodeRightLabelProperty>::const_type
		nodeRightLabel = boost::get( Rule::NodeRightLabelProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::EdgeContextProperty>::const_type
		edgeContext = boost::get( Rule::EdgeContextProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::EdgeLabelProperty>::const_type
		edgeLabel = boost::get( Rule::EdgeLabelProperty(), graph );

	typedef
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map< Rule::CoreGraph::vertex_descriptor, size_t>
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map< Rule::CoreGraph::vertex_descriptor, size_t>
#else
		std::map< Rule::CoreGraph::vertex_descriptor, size_t>
#endif
		Node2IndexMap;

	Node2IndexMap node2idx;

	// write nodes
	Rule::CoreGraph::vertex_iterator vi, v_end;
	boost::tie(vi,v_end) = boost::vertices(graph);
	size_t idx = 0;
	while (vi != v_end) {
		node2idx[*vi] = idx;
		// node head
		out	<<(withSpaces?"  ":"")
			<<"node"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?" ":"")
			<<"id "
			<<idx
			<<" label \""
			;
		// get according label
		switch (nodeContext[*vi]) {
		case Rule::RULE_LEFT_SIDE :
				out <<nodeLabel[*vi] <<"|";
			break;
		case Rule::RULE_RIGHT_SIDE :
				out <<"|" <<nodeLabel[*vi];
			break;
		case Rule::RULE_LABEL_CHANGE :
				out <<nodeLabel[*vi] <<"|" <<nodeRightLabel[*vi];
			break;
		case Rule::RULE_CONTEXT :
				out <<nodeLabel[*vi];
			break;
		}
		// node tail
		out <<"\""
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;

		idx++;
		vi++;
	}


	typedef
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map< std::string, std::string >
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map< std::string, std::string >
#else
		std::map< std::string, std::string >
#endif
		EdgeLabelMap;

	EdgeLabelMap edgeLeft, edgeRight;

	// get edgedata
	Rule::CoreGraph::edge_iterator ei, e_end;
	for (boost::tie(ei,e_end) = boost::edges(graph); ei != e_end; ++ei) {

		std::stringstream edgeID;
		if (node2idx[boost::source(*ei,graph)] <= node2idx[boost::target(*ei,graph)]) {
			edgeID <<(withSpaces?" ":"")
					<<"source "
					<<node2idx[boost::source(*ei,graph)]
					<<" target "
					<<node2idx[boost::target(*ei,graph)];
		} else {
			edgeID <<(withSpaces?" ":"")
							<<"source "
							<<node2idx[boost::target(*ei,graph)]
							<<" target "
							<<node2idx[boost::source(*ei,graph)];
		}

		// get according label
		switch (edgeContext[*ei]) {
		case Rule::RULE_LEFT_SIDE :
				edgeLeft[edgeID.str()] = edgeLabel[*ei];
			break;
		case Rule::RULE_RIGHT_SIDE :
				edgeRight[edgeID.str()] = edgeLabel[*ei];
			break;
		case Rule::RULE_CONTEXT :
				edgeLeft[edgeID.str()] = edgeLabel[*ei];
				edgeRight[edgeID.str()] = edgeLabel[*ei];
			break;
		default :
			throw std::runtime_error("Rule_GML_writer.writeCompact() : edge type 'label change' unexpected not handled");
		}

	}

	// handle left, context, and label change edges
	for (EdgeLabelMap::const_iterator cEdge = edgeLeft.begin(); cEdge!=edgeLeft.end(); ++cEdge) {
		out	<<(withSpaces?"  ":"")
			<<"edge"
			<<(withSpaces?" ":"")
			<<"["
			<< cEdge->first
			<<" label \"";
		if (edgeRight.find(cEdge->first) == edgeRight.end()) {
			// left side edge
			out	<<cEdge->second <<"|";
		} else {
			EdgeLabelMap::const_iterator rightLabel = edgeRight.find(cEdge->first);
			if ( rightLabel->second == cEdge->second ) {
				// context since same label left and right
				out	<<cEdge->second;
			} else {
				// label change
				out <<cEdge->second <<"|" <<rightLabel->second;
			}
			// remove from right side list
			edgeRight.erase( rightLabel );
		}
		out <<"\""
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
		ei++;
	}

	// handle right edges
	for (EdgeLabelMap::const_iterator cEdge = edgeRight.begin(); cEdge!=edgeRight.end(); ++cEdge) {
		out	<<(withSpaces?"  ":"")
			<<"edge"
			<<(withSpaces?" ":"")
			<<"["
			<<(withSpaces?" ":"")
			<< cEdge->first
			<<" label \"";
			// only right side only edges left
		out	<<"|" <<cEdge->second;
		out <<"\""
			<<(withSpaces?" ":"")
			<<"]"
			<<(withSpaces?"\n":"")
			;
		ei++;
	}


	 // close GML output
	out <<"]"
		<<(withSpaces?"\n":"")
		;

}

////////////////////////////////////////////////////////////////////////////////

std::string
Rule_GML_writer::
getContextGML( const Rule::CoreGraph & graph
			, const Rule::RuleContext & context
			, const bool withSpaces )
{

	boost::property_map<Rule::CoreGraph, Rule::NodeContextProperty>::const_type
		nodeContext = boost::get( Rule::NodeContextProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::NodeLabelProperty>::const_type
		nodeLabel = boost::get( Rule::NodeLabelProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::NodeRightLabelProperty>::const_type
		nodeRightLabel = boost::get( Rule::NodeRightLabelProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::EdgeContextProperty>::const_type
		edgeContext = boost::get( Rule::EdgeContextProperty(), graph );
	boost::property_map<Rule::CoreGraph, Rule::EdgeLabelProperty>::const_type
		edgeLabel = boost::get( Rule::EdgeLabelProperty(), graph );

	typedef
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map< Rule::CoreGraph::vertex_descriptor, size_t>
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map< Rule::CoreGraph::vertex_descriptor, size_t>
#else
		std::map< Rule::CoreGraph::vertex_descriptor, size_t>
#endif
		Node2IndexMap;

	Node2IndexMap node2idx;

	std::stringstream out;

	  // write nodes
	Rule::CoreGraph::vertex_iterator vi, v_end;
	boost::tie(vi,v_end) = boost::vertices(graph);
	size_t idx = 0;
	while (vi != v_end) {
		node2idx[*vi] = idx;
		
		if ( nodeContext[*vi] == context
			|| (context < Rule::RULE_LABEL_CHANGE && nodeContext[*vi]==Rule::RULE_LABEL_CHANGE))
		{
			out	<<(withSpaces?"  ":"")
				<<"node"
				<<(withSpaces?" ":"")
				<<"["
				<<(withSpaces?" ":"")
				<<"id "
				<<idx
				<<" label \"";
			  // print according label
			if (nodeContext[*vi] == context || context == Rule::RULE_LEFT_SIDE) {
				out <<nodeLabel[*vi];
			} else {
				out <<nodeRightLabel[*vi];
			}
			out <<"\""
				<<(withSpaces?" ":"")
				<<"]"
				<<(withSpaces?"\n":"")
				;
		}

		idx++;
		vi++;
	}
	
	// write edges
	Rule::CoreGraph::edge_iterator ei, e_end;
	boost::tie(ei,e_end) = boost::edges(graph);
	while (ei != e_end) {
		
		if ( edgeContext[*ei] == context ) {
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
		}
		ei++;
	}

	return out.str();
}


} // namespace ggl
