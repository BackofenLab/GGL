

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

#include <boost/graph/copy.hpp>

#include <iterator>
#include <algorithm>
#include <cassert>

#include <ggl/MR_ApplyRule.hh>

namespace ggl {

////////////////////////////////////////////////////////////////////////////////


	void
	MR_ApplyRule
	::reportHit (	const sgm::Pattern_Interface & pattern,
					const sgm::Graph_Interface & target,
					const sgm::Match & match )
	{
		
#if HAVE_UNORDERED_MAP > 0
//		typedef std::unordered_map<int,sgm::Graph_Interface::IndexType> Label2ShiftMap;
		typedef std::unordered_map<size_t, size_t> StStMap;
#elif HAVE_TR1_UNORDERED_MAP > 0
//		typedef std::tr1::unordered_map<int,sgm::Graph_Interface::IndexType> Label2ShiftMap;
		typedef std::tr1::unordered_map<size_t, size_t> StStMap;
#elif HAVE_GNU_HASH_MAP > 0
		class hash_int {
		public:

			size_t operator()(const int& v) const
			{
				return (size_t)v;
			}
			     
		};
//		typedef __gnu_cxx::hash_map<int,sgm::Graph_Interface::IndexType,sgm::hash_int> Label2ShiftMap;
		typedef __gnu_cxx::hash_map<size_t, size_t> StStMap;
#else
//		typedef std::map<int,sgm::Graph_Interface::IndexType> Label2ShiftMap;
		typedef std::map<size_t, size_t> StStMap;
#endif
		
//		typedef StStMap Graph2SizeMap;

		
		  // check for cast success
		assert( dynamic_cast<const LeftSidePattern*>(&pattern) != NULL /* pattern is no LeftSidePattern object ! */);
		  // cast graph to GGL rule
		const LeftSidePattern* leftSide
			= static_cast<const LeftSidePattern*>(&pattern);


		////////////// PREPARE RESULT GRAPH TO BE FILLED /////////////////////

		  // get first index of each component (direct access)
		typedef ggl::LeftSidePattern::IndexSet IdxSet;
		const IdxSet & firstOfEach = leftSide->getFirstOfEachComponent();
		
		
		  // access to the rule that was matched
		const Rule & rule = leftSide->getRule();
		
		  // check if the left side was really matched
		assert( rule.getLeftSide().nodes.size() == match.size() /*left side and match differ in size*/);
		
		  // the graph that will contain the result of the rule application
		Graph result;

		  // property access maps for the result graph
		boost::property_map<	Graph, PropNodeLabel >
			::type resultNodeLabel = boost::get( 
								PropNodeLabel(), result );
 
		boost::property_map<	Graph, PropNodeIndex >
			::type resultNodeIndex = boost::get( 
								PropNodeIndex(), result );
 
		boost::property_map<	Graph, PropEdgeLabel >
			::type resultEdgeLabel = boost::get( 
								PropEdgeLabel(), result );
 		
		
		  // the mapping of molecule node indices to ITS node indices
		StStMap t2r;
		typedef StStMap::value_type M2Ival;
		

		/////////////  ADD PATTERN NODES (RULE LEFTSIDE) TO RESULT  ///////
		
		  // add all nodes from left side pattern
		for (size_t n = 0; n<leftSide->getNodeNumber(); ++n) {
			  // add new vertex
			Graph::vertex_descriptor newNode = boost::add_vertex(result);
			  // set label
			resultNodeLabel[newNode] = target.getNodeLabel(match.at(n));
			  // store mapping
			t2r.insert(M2Ival(match.at(n),resultNodeIndex[newNode]));
		}


		  // add ALL EDGES between pattern matched nodes, independently if
		  // part of different components im pattern
		if (!addEachComponent) {
			  // check if each node is mapped onto a unique one in target
			if ( t2r.size() != leftSide->getNodeNumber() ) {
				/*************************************************************
				 * some nodes of the leftside rule are mapped onto the same
				 * node in the target graph and no separate copies have to
				 * be created
				 * --> aborting rule application for this match
				 ************************************************************/
				return;
			}
			  // container for special handling of loop labels
			std::multiset<std::string> loopLabels;
			  // add all edges using the target graph
			for (size_t n = 0; n<match.size(); ++n) {
				  // clear label container for next iteration
				loopLabels.clear();
				  // the mapped index of the beginning of the edge
				StStMap::const_iterator fromIdx = t2r.find(match.at(n));
				 // iterate through all out edges
				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(match.at(n)),
						curEdgeEnd = target.getOutEdgesEnd(match.at(n));
					curEdge != curEdgeEnd; ++curEdge )
				{
					  // check if from >= to to avoid duplicated insertions
					if ( undirectedRule && curEdge->getToIndex() < curEdge->getFromIndex() )
						continue;
					  // the mapped index of the end of the edge
					StStMap::const_iterator toIdx = t2r.find(curEdge->getToIndex());
					  // check if edge end is among the present nodes 
					if ( toIdx != t2r.end() )
					{
						  // special handling of loops
						if ( fromIdx->second == toIdx->second )
						{
							if (loopLabels.find( curEdge->getEdgeLabel() ) == loopLabels.end()) {
								  // remember that we have seen this loop
								loopLabels.insert( curEdge->getEdgeLabel() );
							} else {
								  // loop label already seen once
								  // --> since each loop is reported twice, we will ignore this out-edge
								loopLabels.erase( loopLabels.find( curEdge->getEdgeLabel() ) );
								continue;
							}
						}
						  // add new edge
						Graph::edge_descriptor newEdge
						= boost::add_edge(	boost::vertex( fromIdx->second, result),
											boost::vertex( toIdx->second,   result),
											result ).first;
						  // set label
						resultEdgeLabel[newEdge] = curEdge->getEdgeLabel();
					}
				}
			}
		}
		
			
		/////////////  ADD SUB GRAPH(S) TO RESULT  ////////////////////
		
		const sgm::Graph_Interface::CompLabel & compLabel
				= leftSide->getComponentLabeling();

		  // add the matched graph of each component to result
		for (	IdxSet::const_iterator it = firstOfEach.begin();
				it!=firstOfEach.end(); ++it ) 
		{
			  // if each component should have an own molecule copy, 
			  // add node indices of this component to the mapping
			if (addEachComponent) {
				const size_t curLabel = compLabel.at(*it);
				t2r.clear();
				size_t nodesInComponent = 0;
				  // add all index matches of THIS component
				for (size_t i=0; i<match.size(); ++i ) {
					if (compLabel.at(i) == curLabel) {
						assert( t2r.find(match.at(i))==t2r.end() /* two left side nodes of same component are matched onto one target node*/);
						t2r.insert(M2Ival(match.at(i),i));
						nodesInComponent++;
					}
				}
				  // check if each node is mapped onto a unique one in target
				if ( t2r.size() != nodesInComponent ) {
					/*************************************************************
					 * some nodes of the leftside rule are mapped onto the same
					 * node in the target graph and no separate copies have to
					 * be created
					 * --> aborting rule application for this match
					 ************************************************************/
					return;
				}

				  // container for special handling of loop labels
				std::multiset< std::string > loopLabels;
				  // add all edges of this component using the target graph
				for (StStMap::const_iterator fromIdx = t2r.begin(); fromIdx!=t2r.end(); ++fromIdx) {
					  // clear container for next iteration
					loopLabels.clear();
					 // iterate through all adjacent nodes
					for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(fromIdx->first),
							curEdgeEnd = target.getOutEdgesEnd(fromIdx->first);
						curEdge != curEdgeEnd; ++curEdge )
					{
						  // check if index is higher to avoid duplicated insertions
						if ( undirectedRule && curEdge->getToIndex() < curEdge->getFromIndex())
							continue;
						  // the mapped index of the end of the edge
						StStMap::const_iterator toIdx = t2r.find(curEdge->getToIndex());
						  // check if edge end is among the present nodes 
						if ( toIdx != t2r.end() )
						{
							  // special handling of loops
							if ( fromIdx->second == toIdx->second )
							{
								if (loopLabels.find( curEdge->getEdgeLabel() ) == loopLabels.end()) {
									  // remember that we have seen this loop
									loopLabels.insert( curEdge->getEdgeLabel() );
								} else {
									  // loop label already seen once
									  // --> since each loop is reported twice, we will ignore this out-edge
									loopLabels.erase( loopLabels.find( curEdge->getEdgeLabel() ) );
									continue;
								}
							}
							  // add new edge
							Graph::edge_descriptor newEdge
							= boost::add_edge(	boost::vertex( fromIdx->second, result),
												boost::vertex( toIdx->second,   result),
												result ).first;
							  // set label
							resultEdgeLabel[newEdge] = curEdge->getEdgeLabel();
						}
					}
				}
				
			} // if addEachComponent
			
			  // will contain all indices of target that are part of the
			  // molecule component matched by the current rule component
			std::set< size_t > nodesToHandle, handledNodes;
			  // add first of this component
			nodesToHandle.insert(match.at(*it));
			
			  // traverse current molecule graph recursively from node to node
			while (!nodesToHandle.empty()) {
				  // get first non-handled node index and remove from storage
				const size_t curNodeIdx = *(nodesToHandle.begin());
				nodesToHandle.erase(nodesToHandle.begin());
				  // remember node as already handled
				handledNodes.insert(curNodeIdx);
				
				const bool isNewNode = t2r.find(curNodeIdx) == t2r.end();
				
				  // check if not already part of the result graph
				if (isNewNode) {
					  // add new vertex
					Graph::vertex_descriptor newNode = boost::add_vertex(result);
					  // store mapping
					t2r.insert(M2Ival( curNodeIdx, resultNodeIndex[newNode]));
					  // set label
					resultNodeLabel[newNode] = target.getNodeLabel(curNodeIdx);
				}
				
				  // container for special handling of loops
				std::multiset< std::string > loopLabels;
				 // iterate through all adjacent nodes
				for( sgm::Graph_Interface::OutEdge_iterator curEdge = target.getOutEdgesBegin(curNodeIdx),
						curEdgeEnd = target.getOutEdgesEnd(curNodeIdx);
					curEdge != curEdgeEnd; ++curEdge )
				{
					
					assert( curEdge->getFromIndex() == curNodeIdx /*otherwise the following code is incomplete*/);
					
					  // if neighbored node was not already handled add to
					  // processing queue
					if (handledNodes.find(curEdge->getToIndex())==handledNodes.end()) {
						nodesToHandle.insert(curEdge->getToIndex());
					}

					StStMap::const_iterator fromIdx = t2r.find(curEdge->getFromIndex());
					StStMap::const_iterator toIdx = t2r.find(curEdge->getToIndex());
					
					  // add edge if at least one new node involved and both already added
					if (	toIdx != t2r.end() // == both nodes already added
							&& isNewNode // == new node from target 
						) 
					{
						  // special handling of loops
						if ( fromIdx->second == toIdx->second )
						{
							if (loopLabels.find( curEdge->getEdgeLabel() ) == loopLabels.end()) {
								  // remember that we have seen this loop
								loopLabels.insert( curEdge->getEdgeLabel() );
							} else {
								  // loop label already seen once
								  // --> since each loop is reported twice, we will ignore this out-edge
								loopLabels.erase( loopLabels.find( curEdge->getEdgeLabel() ) );
								continue;
							}
						}
						  // add new edge
						std::pair< Graph::edge_descriptor, bool> newEdge
						= boost::add_edge(	boost::vertex( fromIdx->second, result),
											boost::vertex( toIdx->second,   result),
											result);
						assert( newEdge.second /* otherwise edge to insert was already present */);
						  // set label
						resultEdgeLabel[newEdge.first] = curEdge->getEdgeLabel();
					}
				}
				
			} // while not all node handled
			
		} // for each component
			
		////////////////////////////////////////////////////////////////////////
		// NOW ALL NODES AND EDGES OF THE MAPPED GRAPH (COMPONENT) ARE COPIED
		////////////////////////////////////////////////////////////////////////

		StStMap::const_iterator t2r_const_end = t2r.end();
			
		////////////////// APPLY RULE ///////////////////////
		
		  // get access to the nodes
		typedef Rule::CoreGraph::vertex_descriptor NodeDescr;
//		typedef Rule::CoreGraph::edge_descriptor EdgeDescr;

		  // get access to the necessary property maps of the rule
		
		boost::property_map<	Rule::CoreGraph
								, Rule::EdgeLabelProperty
							>
			::const_type ruleEdgeLabel = boost::get(	
								Rule::EdgeLabelProperty()
								, rule.getCore() );

		boost::property_map<	Rule::CoreGraph
								, Rule::NodeLabelProperty
							>
			::const_type ruleNodeLabel = boost::get(	
								Rule::NodeLabelProperty()
								, rule.getCore() );
		boost::property_map<	Rule::CoreGraph
								, Rule::NodeRightLabelProperty
							>
			::const_type ruleNodeRightLabel = boost::get(	
								Rule::NodeRightLabelProperty()
								, rule.getCore() );
		boost::property_map<	Rule::CoreGraph
								, Rule::NodeIndexProperty
							>
			::const_type ruleNodeIndex = boost::get(	
								Rule::NodeIndexProperty()
								, rule.getCore() );

		boost::property_map<	Rule::CoreGraph
								, Rule::NodeContextProperty
							>
			::const_type ruleNodeContext = boost::get(	
								Rule::NodeContextProperty()
								, rule.getCore() );

		boost::property_map<	Rule::CoreGraph
								, Rule::EdgeContextProperty
							>
			::const_type ruleEdgeContext = boost::get(	
								Rule::EdgeContextProperty()
								, rule.getCore() );
		
		  // set up node index mapping for rule graph
		t2r.clear();
		for (size_t i=0; i<rule.getLeftSide().nodes.size(); ++i ) {
			  // add node mapping for already inserted left side rule nodes
			t2r.insert(M2Ival( rule.getLeftSide().nodes.at(i), i));
		}

		////////////////////////////////////////////////////////////////////////
		// add new and rename nodes
		////////////////////////////////////////////////////////////////////////
		
		 // check for renaming and adding of nodes (need right-side only)
		for (size_t i=0; i<rule.getRightSide().nodes.size(); ++i ) {
			  // the index of the i-th right side node within the CoreGraph
			const size_t ruleNodeIdx = rule.getRightSide().nodes.at(i);
			  // access to the rule node
			NodeDescr node = boost::vertex(	ruleNodeIdx, 
											rule.getCore());
			  // change result graph according to the node context
			switch ( ruleNodeContext[node] ) 
			{
				  // rename target node
				case ggl::Rule::RULE_LABEL_CHANGE : {
					  // get access to matched node 
					assert(t2r.find(ruleNodeIdx)!=t2r_const_end /* otherwise node to change label of doesnt exist */);
					Graph::vertex_descriptor
						rNode = boost::vertex(
								t2r.find( ruleNodeIdx )->second
					            , result);
					  // change label of matched node
					resultNodeLabel[rNode] = getAlteredNodeLabel( resultNodeLabel[rNode], ruleNodeRightLabel[node] );
					break;
				}
				  // add new target node
				case ggl::Rule::RULE_RIGHT_SIDE : {
					  // add new vertex
					Graph::vertex_descriptor
						rNode = boost::add_vertex(result);
					  // add new node to mapping
					t2r.insert(M2Ival( ruleNodeIdx, resultNodeIndex[rNode]));
					  // set label of new node
					resultNodeLabel[rNode] = ruleNodeLabel[node];
					break;
				}
				default : break;
			}
		}

		////////////////////////////////////////////////////////////////////////
		// NOW ALL NODES OF THE RULE ARE ADDED AND MAPPED TO NODES IN THE GRAPH
		////////////////////////////////////////////////////////////////////////


		  // iterators for graph descriptor access
		Rule::CoreGraph::out_edge_iterator ei, e_end;
		Rule::CoreGraph::vertex_iterator ni, n_end;

		////////////////////////////////////////////////////////////////////////
		// perform copy-and-paste operations
		////////////////////////////////////////////////////////////////////////


		  // apply copy-and-paste operations
		for (Rule::CopyAndPasteOperations::const_iterator op = rule.getCopyAndPasteOperations().begin();
				op != rule.getCopyAndPasteOperations().end(); ++ op)
		{
			assert(t2r.find( op->first ) != t2r_const_end); // otherwise source node unknown
			  // access to the source node
			Graph::vertex_descriptor sourceNode
				= boost::vertex( t2r.find( op->first )->second , result);
			  // perform all paste operations for this source node
			  // get access to the list
			typedef Rule::CopyAndPasteOperations::mapped_type CNPlist;
			const CNPlist & cnpList = op->second;
			  // container for special handling of loop labels
			std::multiset<std::string> loopLabels;
			  // process CnP list
			for (CNPlist::const_iterator it=cnpList.begin(); it!=cnpList.end(); ++it) {
				  // access to the current CnP operation to perform
				const Rule::RuleCnP & curCnP = *it;
				  // assure that the source node is correct
				assert( sourceNode == boost::vertex( t2r.find(curCnP.source)->second, result) );
				assert( t2r.find(curCnP.pasteID) != t2r_const_end );
				  // get paste vertex
				Graph::vertex_descriptor node2paste =
						boost::vertex( t2r.find(curCnP.pasteID)->second, result);
				  // copy all edge labels to cnp pasteID
				  // (no edge labels given
				  //  or the wildcard string is among the edge labels)
				const bool copyAllLabels = (curCnP.edgeLabels.empty()
							|| (pattern.getUsedWildcard()!=NULL
								&& curCnP.edgeLabels.find(*(pattern.getUsedWildcard())) != curCnP.edgeLabels.end()));
				Graph::out_edge_iterator ei,e_end;
				  // clear container for next iteration
				loopLabels.clear();
				  // check all edges
				for (boost::tie(ei,e_end) = out_edges( sourceNode, result );
						ei != e_end; ++ei)
				{

					  // edge label check
					bool copyThisEdge = copyAllLabels ||
							curCnP.edgeLabels.find(resultEdgeLabel[*ei]) != curCnP.edgeLabels.end();

					  // check target node id if necessary
					if (copyThisEdge && curCnP.target != (size_t)INT_MAX ) {
						assert( t2r.find(curCnP.target) != t2r_const_end );
						  // check if target id correct
						copyThisEdge = resultNodeIndex[boost::target(*ei,result)] == t2r.find(curCnP.target)->second;
					}

					  // add according edge if to be copied
					if ( copyThisEdge ) {
						  // special handling of loops
						if ( sourceNode == boost::target(*ei,result) )
						{
							if (loopLabels.find( resultEdgeLabel[*ei] ) == loopLabels.end()) {
								  // remember that we have seen this loop
								loopLabels.insert( resultEdgeLabel[*ei] );
								  // add loop and set label
								Graph::edge_descriptor newEdge
									= boost::add_edge( node2paste, node2paste, result).first;
								resultEdgeLabel[ newEdge ] = resultEdgeLabel[ *ei ];
							} else {
								  // loop label already seen once
								  // --> since each loop is reported twice, we will ignore this out-edge
								loopLabels.erase( loopLabels.find( resultEdgeLabel[*ei] ) );
								continue;
							}
						} else {
							  // add edge and set label
							Graph::edge_descriptor newEdge
								= boost::add_edge( node2paste, boost::target(*ei,result), result).first;
							resultEdgeLabel[ newEdge ] = resultEdgeLabel[ *ei ];
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////
		// add new, remove old, and rename edges
		////////////////////////////////////////////////////////////////////////

		  // access to the wildcard used
		const std::string * wildcard = pattern.getUsedWildcard();

		  // iterate over all nodes and their adjacent outedges
		for ( boost::tie(ni, n_end) = boost::vertices(rule.getCore()); ni != n_end; ++ni ) {

			  // count how many edges per target have to be removed based on wildcard mapping
			std::vector< size_t > removeWildcardEdges(boost::num_vertices(result),0);
			typedef std::multiset< std::string > EdgeLabelSet;
			typedef std::vector< EdgeLabelSet * > EdgeLabelSetVec;
			EdgeLabelSetVec contextEdges(boost::num_vertices(result),NULL);

			assert( t2r.find( ruleNodeIndex[*ni] ) != t2r_const_end /* source node to connect is missing*/);
			size_t rFrom = t2r.find( ruleNodeIndex[*ni] )->second;

			  // check all out edges if they have to be removed
			for (boost::tie(ei, e_end) = boost::out_edges(*ni, rule.getCore()); ei != e_end; ++ei) {

				  // get edge node indices of the edge in the result/target graph
				assert( t2r.find( ruleNodeIndex[boost::target(*ei,rule.getCore())] )!=t2r_const_end /* target node to connect is missing*/);
				size_t rTo   = t2r.find( ruleNodeIndex[boost::target(*ei,rule.getCore())] )->second;
				  // order check to avoid double processing
				if ( undirectedRule && rTo < rFrom ) {
					continue;
				}

				  // context edge handling
				if ( wildcard != NULL && ruleEdgeContext[*ei] == ggl::Rule::RULE_CONTEXT
						&& ruleEdgeLabel[*ei].compare(*wildcard) != 0 ) // store only non-wildcard context label
				{
					  // store context edge label to ensure that it is not removed by wildcard matches
					if (contextEdges[rTo] == NULL)
						contextEdges[rTo] = new EdgeLabelSet();
					contextEdges[rTo]->insert(ruleEdgeLabel[*ei]);
					continue;
				}
				  // the remaining handling is for left side edges only
				if ( ruleEdgeContext[*ei] != ggl::Rule::RULE_LEFT_SIDE ) {
					continue;
				}
				  // wildcard handling
				if ( wildcard != NULL && ruleEdgeLabel[*ei].compare(*wildcard)==0 ) {
					  // count wildcard removal
					removeWildcardEdges[ rTo ] += 1;
				} else {
					  // find edge with according label in result graph and remove
					Graph::out_edge_iterator tei, te_end;
					bool edgeNotFound = true;
					for (boost::tie(tei, te_end) = boost::out_edges(boost::vertex(rFrom,result), result); edgeNotFound && tei != te_end; ++tei) {
						  // check if correct edge found
						if ( resultNodeIndex[boost::target( *tei, result )] == rTo
							&& resultEdgeLabel[*tei].compare(ruleEdgeLabel[*ei]) == 0 )
						{
							  // remove edge
							edgeNotFound = false;
							boost::remove_edge( *tei, result );
						}
					}
					assert( !edgeNotFound /* otherwise edge was not found */);
				}

			} // all outedges of current node

			// remove wildcard matches but take care about context stuff that has to be maintained
			if (wildcard != NULL) {
				  // get overall number to be removed
				size_t edgesToRemove = 0;
				for (size_t i=0; i<removeWildcardEdges.size(); ++i) {
					edgesToRemove += removeWildcardEdges.at(i);
				}
				if (edgesToRemove > 0) {
					  // collect candidate edges for removal
					Graph::out_edge_iterator tei, te_end;
					EdgeLabelSetVec currentLabel(boost::num_vertices(result),NULL);
					for (boost::tie(tei, te_end) = boost::out_edges(boost::vertex(rFrom,result), result); tei != te_end; ++tei) {
						  // check if we have to remove an edge to the current target
						size_t rTo = resultNodeIndex[ boost::target(*tei,result) ];
						if ( removeWildcardEdges.at(rTo) > 0 ) {
							 // store edge label for later comparison
							if (currentLabel.at(rTo) == NULL)
								currentLabel[rTo] = new EdgeLabelSet();
							currentLabel[rTo]->insert( resultEdgeLabel[*tei] );
						}
					}
					  // identify labels that can be removed
					for (size_t i=0; i<currentLabel.size(); ++i) {
						   // check if context label available
						if (currentLabel.at(i) != NULL && contextEdges.at(i) != NULL) {
							  // remove context label from current labels
							std::multiset<std::string> diff;
							std::insert_iterator< std::multiset<std::string> > diffInserter(diff,diff.begin());
							std::set_difference( currentLabel.at(i)->begin(), currentLabel.at(i)->end(),
													contextEdges.at(i)->begin(), contextEdges.at(i)->end(),
													diffInserter );
							assert(diff.size() >= removeWildcardEdges.at(i) /* otherwise less labels available than to be removed */);
							  // store labels we are allowed to remove
							*(currentLabel[i]) = diff;
						}
					}
					bool edgeRemoved = false;
					do {
						edgeRemoved = false;
						for (boost::tie(tei, te_end) = boost::out_edges(boost::vertex(rFrom,result), result); !edgeRemoved && tei != te_end; ++tei) {
							  // check if we have to remove an edge to the current target
							size_t rTo = resultNodeIndex[ boost::target(*tei,result) ];
							if ( removeWildcardEdges.at(rTo) > 0 ) {
								assert( currentLabel.at(rTo)!=NULL);
								if( currentLabel.at(rTo)->find( resultEdgeLabel[*tei] ) != currentLabel.at(rTo)->end() )
								{
									  // remove edge information
									currentLabel[rTo]->erase(currentLabel.at(rTo)->find( resultEdgeLabel[*tei] ));
									  // decrease counter
									removeWildcardEdges[rTo] -= 1;
									--edgesToRemove;
									  // special handling of loops
									  // --> if you remove one loop you remove two out edges
									if (rFrom == rTo) {
										  // remove back edge information
										currentLabel[rTo]->erase(currentLabel.at(rTo)->find( resultEdgeLabel[*tei] ));
										  // decrease counter again
										removeWildcardEdges[rTo] -= 1;
										--edgesToRemove;
									}
									  // remove edge
									boost::remove_edge( *tei, result );
									  // restart iteration since iterators might be invalidated
									edgeRemoved = true;
								}
							}
						}
					} while( edgeRemoved );
					assert( edgesToRemove == 0 /*somehow we could not delete enough edges*/ );

					  // clear temporary data
					for (EdgeLabelSetVec::iterator l=currentLabel.begin(); l!=currentLabel.end(); ++l) {
						if (*l != NULL)
							delete *l;
					}
				}

			}

			  // container for special handling of loops
			std::multiset<std::string> loopLabels;
			// add new edges from right side context
			for (boost::tie(ei, e_end) = boost::out_edges(*ni, rule.getCore()); ei != e_end; ++ei) {
				  // skip other edges in case no wildcard used
				if ( ruleEdgeContext[*ei] != ggl::Rule::RULE_RIGHT_SIDE ) {
					continue;
				}

				  // get edge node indices of the edge in the result/target graph
				assert( t2r.find( ruleNodeIndex[boost::source(*ei,rule.getCore())] )!=t2r_const_end /* source node to connect is missing*/);
				assert( t2r.find( ruleNodeIndex[boost::target(*ei,rule.getCore())] )!=t2r_const_end /* target node to connect is missing*/);
				size_t rFrom = t2r.find( ruleNodeIndex[boost::source(*ei,rule.getCore())] )->second;
				size_t rTo   = t2r.find( ruleNodeIndex[boost::target(*ei,rule.getCore())] )->second;
				  // order check to avoid double processing
				if ( undirectedRule && rTo < rFrom ) {
					continue;
				}
				  // special handling of loops
				if ( rFrom == rTo )
				{
					if (loopLabels.find( ruleEdgeLabel[*ei] ) == loopLabels.end()) {
						  // remember that we have seen this loop
						loopLabels.insert( ruleEdgeLabel[*ei] );
					} else {
						  // loop label already seen once
						  // --> since each loop is reported twice, we will ignore this out-edge
						loopLabels.erase( loopLabels.find( ruleEdgeLabel[*ei] ) );
						continue;
					}
				}

				  // add new edge
				std::pair<Graph::edge_descriptor, bool> rEdge
					= boost::add_edge(	boost::vertex( rFrom, result),
										boost::vertex( rTo,   result),
										result);
				  // set label
				resultEdgeLabel[rEdge.first] = ruleEdgeLabel[*ei];
			}

			  // clear temporary data
			for (EdgeLabelSetVec::iterator l=contextEdges.begin(); l!=contextEdges.end(); ++l) {
				if (*l != NULL)
					delete *l;
			}

		} // all nodes

		////////////////////////////////////////////////////////////////////////
		// remove nodes
		////////////////////////////////////////////////////////////////////////

		  // a list of all node indices in result that have to be removed
		  // --> use original rule indices !!!
		std::vector<size_t> toRemove;
		
		 // check for nodes to remove (need left-side only)
		for (size_t i=0; i<rule.getLeftSide().nodes.size(); ++i ) {
			  // access to the rule node
			NodeDescr node = boost::vertex(	rule.getLeftSide().nodes.at(i), 
											rule.getCore());
			  // change result graph according to the node context
			if ( ruleNodeContext[node] == ggl::Rule::RULE_LEFT_SIDE ) {
				  // store index to remove 
				assert(t2r.find(rule.getLeftSide().nodes.at(i)) != t2r_const_end);
				toRemove.push_back( t2r.find(rule.getLeftSide().nodes.at(i))->second );
			}
		} // now we know all nodes that have to be removed
		
		  // sort the indices to remove
		std::sort(toRemove.begin(), toRemove.end());
		
		  // remove nodes in decreasing index order
		  //  to avoid iterator invalidation
		for (std::vector<size_t>::const_reverse_iterator remFrom =toRemove.rbegin();
				remFrom != toRemove.rend(); ++remFrom)
		{
			  // access node to delete
			Graph::vertex_descriptor node2delete
				= boost::vertex( *remFrom , result);
			while( boost::edge(node2delete, node2delete, result).second) {
			}
			  // remove adjacent edges
			boost::clear_vertex( node2delete, result);
			  // check if loop present (bug in boost < 1.50.0)
			if( boost::edge(node2delete, node2delete, result).second) {
				boost::remove_edge(node2delete,node2delete,result);
			}
			  // remove the node
			boost::remove_vertex( node2delete, result);
		}


		  // add result to storage
		storage.add(result);
		
	}
	
////////////////////////////////////////////////////////////////////////////////

	
	
////////////////////////////////////////////////////////////////////////////////

} // namespace ggl

