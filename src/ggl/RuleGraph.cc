
#include <cassert>

#include <sgm/Graph_Interface.hh>
#include <ggl/RuleGraph.hh>

namespace ggl {

////////////////////////////////////////////////////////////////////////////////


	LeftSidePattern
	::LeftSidePattern (const Rule& rule_)
	 :	rule(rule_)
	 	, firstOfEachComponent()
	 	, compLabels()
	 	, symmBreakList(NULL)
		, edgeLabel(boost::get( Rule::EdgeLabelProperty(), rule.getCore() ))
		, nodeLabel(boost::get( Rule::NodeLabelProperty(), rule.getCore() ))
		, nodeIndex(boost::get( Rule::NodeIndexProperty(), rule.getCore() ))
		, edgeContext(boost::get( Rule::EdgeContextProperty(), rule.getCore() ))
		, nodeContext(boost::get( Rule::NodeContextProperty(), rule.getCore() ))
	{
		// get labeling of components
		sgm::Graph_Interface::connectedComponents( *this, compLabels );
		// get first index of each component
		std::set<int> doneComponents;
		for (size_t i=0; i<compLabels.size(); ++i ) {
			  // check if component was already seen
			if ( doneComponents.find(compLabels[i]) == doneComponents.end() ) {
				  // new component --> store first
				firstOfEachComponent.insert(i);
				  // store that this component label was seen
				doneComponents.insert(compLabels[i]);
			}
		}

	}


////////////////////////////////////////////////////////////////////////////////


	LeftSidePattern
	::LeftSidePattern (const LeftSidePattern& toCopy )
	 :	rule( toCopy.rule )
	 	, firstOfEachComponent( toCopy.firstOfEachComponent )
	 	, compLabels( toCopy.compLabels )
	 	, symmBreakList( toCopy.symmBreakList == NULL
 					? NULL
 					: new sgm::PA_OrderCheck::CheckList(*(toCopy.symmBreakList)) )
		, edgeLabel(boost::get( Rule::EdgeLabelProperty(), rule.getCore() ))
		, nodeLabel(boost::get( Rule::NodeLabelProperty(), rule.getCore() ))
		, nodeIndex(boost::get( Rule::NodeIndexProperty(), rule.getCore() ))
		, edgeContext(boost::get( Rule::EdgeContextProperty(), rule.getCore() ))
		, nodeContext(boost::get( Rule::NodeContextProperty(), rule.getCore() ))
	{
	}


//////////////////////////////////////////////////////////////////////////////////


	sgm::Graph_Interface::OutEdge_iterator
	LeftSidePattern
	::getOutEdgesBegin(const IndexType & i) const
	{
		assert( i<getNodeNumber() );
		  // get edge iterators for given node
		Rule::CoreGraph::out_edge_iterator start,end;
		boost::tie(start,end) = boost::out_edges( boost::vertex(rule.getLeftSide().nodes.at(i),rule.getCore()), rule.getCore());
		  // return according iterator that describes the iteration start
		return OutEdge_iterator( EdgeDescriptor( start, end, *this ) );
	}

////////////////////////////////////////////////////////////////////////////////


	sgm::Graph_Interface::OutEdge_iterator
	LeftSidePattern::
	getOutEdgesEnd( const IndexType & i ) const
	{
		assert( i<getNodeNumber() );
		  // get edge iterators for given node
		Rule::CoreGraph::out_edge_iterator end
			= boost::out_edges( boost::vertex(rule.getLeftSide().nodes.at(i),rule.getCore()), rule.getCore()).second;
		  // return according iterator that only describes the iteration end
		return OutEdge_iterator( EdgeDescriptor( end, end, *this ) );
	}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



	LeftSidePattern::
	EdgeDescriptor::
	EdgeDescriptor(	const Rule::CoreGraph::out_edge_iterator& cur_edge_,
					const Rule::CoreGraph::out_edge_iterator& list_end_,
					const LeftSidePattern & parent_ )
	 : Graph_Interface::EdgeDescriptor()
		, cur_edge( cur_edge_ )
		, list_end( list_end_ )
		, parent(parent_)
	{

		  // skip pure right side edges
		while( cur_edge != list_end
				&& (parent.nodeContext[boost::target(*cur_edge,parent.rule.getCore())] == Rule::RULE_RIGHT_SIDE
						|| parent.edgeContext[*cur_edge] == Rule::RULE_RIGHT_SIDE ))
		{
			++cur_edge;
		}

		  // check if no neighbors available
		if (cur_edge != list_end) {
			assert(parent.nodeContext[boost::source(*cur_edge,parent.rule.getCore())] != Rule::RULE_RIGHT_SIDE /* otherwise right side node*/);
			  // set source node
			from = parent.nodeIndex[ boost::source( *cur_edge, parent.rule.getCore() ) ];
			  // map index to left side pattern graph node index
			from = (size_t)(std::find(	parent.rule.getLeftSide().nodes.begin(),
										parent.rule.getLeftSide().nodes.end(),
										from )
							- parent.rule.getLeftSide().nodes.begin());
			assert( parent.rule.getLeftSide().nodes.at(from) == parent.nodeIndex[ boost::source( *cur_edge, parent.rule.getCore() ) ]);
			  // set target node index within complete rule graph
			to = parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ];
			  // map index to left side pattern graph node index
			to = (size_t)(std::find(	parent.rule.getLeftSide().nodes.begin(),
										parent.rule.getLeftSide().nodes.end(),
										to )
							- parent.rule.getLeftSide().nodes.begin());
			assert( parent.rule.getLeftSide().nodes.at(to) == parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ]);
		}
	}


	//////////////////////////////////////////////////////////////////////////


	LeftSidePattern::
	EdgeDescriptor&
	LeftSidePattern::
	EdgeDescriptor::
	operator++()
	{
		if (cur_edge != list_end) {
			  // increase iterator
			++cur_edge;
			  // skip pure right side edges
			while( cur_edge != list_end
					&& (parent.edgeContext[*cur_edge] == Rule::RULE_RIGHT_SIDE
							|| parent.nodeContext[boost::target(*cur_edge,parent.rule.getCore())] == Rule::RULE_RIGHT_SIDE ))
			{
				++cur_edge;
			}
			  // update data
			  // check if no neighbors available --> set to NO_FURTHER_EDGES
			if (cur_edge != list_end) {
				  // set target node
				to = parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ];
				  // map index to left side pattern graph node index
				to = (size_t)(std::find(	parent.rule.getLeftSide().nodes.begin(),
											parent.rule.getLeftSide().nodes.end(),
											to )
								- parent.rule.getLeftSide().nodes.begin());
				assert( parent.rule.getLeftSide().nodes.at(to) == parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ]);
			}

		}

		return *this;
	}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	
} // namespace ggl


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace ggl {


////////////////////////////////////////////////////////////////////////////////


	std::string
	RightSidePattern
	::getNodeLabel(const IndexType & i) const
	{
		assert( i < rule.getRightSide().nodes.size() /* index out of scope */ );
		
		  // get access to the node i
		Rule::CoreGraph::vertex_descriptor nodeDescr
			= boost::vertex( rule.getRightSide().nodes.at(i), rule.getCore());

		  // get label according to rule behavior
		std::string label = "";
		if ( nodeContext[nodeDescr] == Rule::RULE_LABEL_CHANGE ) {
			label = nodeRightLabel[nodeDescr];
		} else {
			label = nodeLabel[nodeDescr];
		}
		
		return label;
		
	}
	
////////////////////////////////////////////////////////////////////////////////


	sgm::Graph_Interface::OutEdge_iterator
	RightSidePattern
	::getOutEdgesBegin(const IndexType & i) const
	{
		assert( i<getNodeNumber() );
		  // get edge iterators for given node
		Rule::CoreGraph::out_edge_iterator start,end;
		boost::tie(start,end) = boost::out_edges( boost::vertex(rule.getRightSide().nodes.at(i),rule.getCore()), rule.getCore());
		  // return according iterator that describes the iteration start
		return OutEdge_iterator( EdgeDescriptor( start, end, *this ) );
	}

////////////////////////////////////////////////////////////////////////////////


	sgm::Graph_Interface::OutEdge_iterator
	RightSidePattern::
	getOutEdgesEnd( const IndexType & i ) const
	{
		assert( i<getNodeNumber() );
		  // get edge iterators for given node
		Rule::CoreGraph::out_edge_iterator end
			= boost::out_edges( boost::vertex(rule.getRightSide().nodes.at(i),rule.getCore()), rule.getCore()).second;
		  // return according iterator that only describes the iteration end
		return OutEdge_iterator( EdgeDescriptor( end, end, *this ) );
	}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


	RightSidePattern::
	EdgeDescriptor::
	EdgeDescriptor(	const Rule::CoreGraph::out_edge_iterator& cur_edge_,
					const Rule::CoreGraph::out_edge_iterator& list_end_,
					const RightSidePattern & parent_ )
	 : Graph_Interface::EdgeDescriptor()
		, cur_edge( cur_edge_ )
		, list_end( list_end_ )
		, parent(parent_)
	{

		  // skip pure left side edges
		while( cur_edge != list_end
				&& ( parent.edgeContext[*cur_edge] == Rule::RULE_LEFT_SIDE
						|| parent.nodeContext[boost::target(*cur_edge,parent.rule.getCore())] == Rule::RULE_LEFT_SIDE ))
		{
			++cur_edge;
		}

		  // check if no neighbors available
		if (cur_edge != list_end) {
			assert(parent.nodeContext[boost::source(*cur_edge,parent.rule.getCore())] != Rule::RULE_LEFT_SIDE /* otherwise right side node*/);
			  // set source node
			from = parent.nodeIndex[ boost::source( *cur_edge, parent.rule.getCore() ) ];
			  // map index to right side pattern graph node index
			from = (size_t)(std::find(	parent.rule.getRightSide().nodes.begin(),
										parent.rule.getRightSide().nodes.end(),
										from )
							- parent.rule.getRightSide().nodes.begin());
			assert( parent.rule.getRightSide().nodes.at(from) == parent.nodeIndex[ boost::source( *cur_edge, parent.rule.getCore() ) ]);
			  // set target node index within complete rule graph
			to = parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ];
			  // map index to right side pattern graph node index
			to = (size_t)(std::find(	parent.rule.getRightSide().nodes.begin(),
										parent.rule.getRightSide().nodes.end(),
										to )
							- parent.rule.getRightSide().nodes.begin());
			assert( parent.rule.getRightSide().nodes.at(to) == parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ]);
		}
	}

	//////////////////////////////////////////////////////////////////////////


	RightSidePattern::
	EdgeDescriptor&
	RightSidePattern::
	EdgeDescriptor::
	operator++()
	{
		if (cur_edge != list_end) {
			  // increase iterator
			++cur_edge;
			  // skip pure left side edges
			while( cur_edge != list_end
					&& ( parent.edgeContext[*cur_edge] == Rule::RULE_LEFT_SIDE
							|| parent.nodeContext[boost::target(*cur_edge,parent.rule.getCore())] == Rule::RULE_LEFT_SIDE  ))
			{
				++cur_edge;
			}
			  // update data
			  // check if no neighbors available --> set to NO_FURTHER_EDGES
			if (cur_edge != list_end) {
				  // set target node
				to = parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ];
				  // map index to right side pattern graph node index
				to = (size_t)(std::find(	parent.rule.getRightSide().nodes.begin(),
											parent.rule.getRightSide().nodes.end(),
											to )
								- parent.rule.getRightSide().nodes.begin());
				assert( parent.rule.getRightSide().nodes.at(to) == parent.nodeIndex[ boost::target( *cur_edge, parent.rule.getCore() ) ]);
			}

		}

		return *this;
	}

	//////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
	

} // namespace ggl
