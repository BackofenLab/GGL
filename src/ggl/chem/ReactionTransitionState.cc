
#include <cassert>
#include <sstream>
#include <algorithm>
#include <boost/graph/copy.hpp>

#include "ggl/chem/ReactionTransitionState.hh"

namespace ggl {
 namespace chem {


//##############################################################################



	ReactionTransitionState::
	ReactionTransitionState(
						const LeftSidePattern & leftSide
						, const sgm::Graph_Interface & molecules
						, const sgm::Match & match
						, const bool addEachComponent )
 	 : SuperClass( )
	{
		  // create the temporary imaginary transition state object to fill
		ChemRule::CoreGraph its;
		  // merge rule core graph and educts (from molecules) into the its
		if (mergeGraphs( its, leftSide, molecules, match, addEachComponent )) {
			  // create final transition state representation
			SuperClass::initializeTransitionState( its );
		} else {
			  // overwrite with empty graph
			its = ChemRule::CoreGraph();
		}
	}
	

//##############################################################################



	bool
	ReactionTransitionState::
	mergeGraphs(	ChemRule::CoreGraph & its
					, const LeftSidePattern & leftSide
					, const sgm::Graph_Interface & molecules
					, const sgm::Match & match
					, const bool addEachComponent )
	{
 		
 		using namespace ggl;
 		using namespace ggl::chem;

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

		typedef ChemRule::CoreGraph CoreGraph;
		typedef LeftSidePattern::IndexSet IdxSet;

		  // access to the rule that was matched
		const ChemRule & rule = leftSide.getRule();
		
		  // check if the left side was really matched
		assert( rule.getLeftSide().nodes.size() == match.size() /*left side and match differ in size*/);
		
		  // copy the original rule graph into the ITS object
		boost::copy_graph< CoreGraph, CoreGraph >( rule.getCore(), its );

		
		  // get access to the ITS nodes
		typedef CoreGraph::vertex_descriptor NodeDescr;
		typedef CoreGraph::edge_descriptor EdgeDescr;

		  // get access to the necessary property maps of the ITS
		
		boost::property_map<	CoreGraph
								, ChemRule::EdgeLabelProperty
							>
			::type itsEdgeLabel = boost::get(	
								ChemRule::EdgeLabelProperty()
								, its );

		boost::property_map<	CoreGraph
								, ChemRule::NodeLabelProperty
							>
			::type itsNodeLabel = boost::get(	
								ChemRule::NodeLabelProperty()
								, its );
//		boost::property_map<	CoreGraph
//								, ChemRule::NodeRightLabelProperty
//							>
//			::type itsNodeRightLabel = boost::get(
//								ChemRule::NodeRightLabelProperty()
//								, its );
		boost::property_map<	CoreGraph
								, ChemRule::NodeIndexProperty
							>
			::type itsNodeIndex = boost::get(	
								ChemRule::NodeIndexProperty()
								, its );

		boost::property_map<	CoreGraph
								, ChemRule::NodeContextProperty
							>
			::type itsNodeContext = boost::get(	
								ChemRule::NodeContextProperty()
								, its );

		boost::property_map<	CoreGraph
								, ChemRule::EdgeContextProperty
							>
			::type itsEdgeContext = boost::get(	
								ChemRule::EdgeContextProperty()
								, its );

		
		  // the mapping of molecule node indices to ITS node indices
		StStMap m2i;
		typedef StStMap::value_type M2Ival;
		
		  // if molecule is to use as it is ...
		if (!addEachComponent) {
			  // add all matches to mapping, independent from components
			for (size_t i=0; i<match.size(); ++i ) {
				m2i.insert(M2Ival(match.at(i),rule.getLeftSide().nodes.at(i)));
			}
			  // check if each node is mapped onto a unique one in target
			if ( m2i.size() != match.size() ) {
				/*************************************************************
				 * some nodes of the leftside rule are mapped onto the same
				 * node in the target graph and no separate copies have to
				 * be created
				 * --> aborting rule application for this match
				 ************************************************************/
				return false;
			}

			// add missing edges from molecule between nodes from the
			// leftside rule pattern, and each only once
			for (size_t n = 0; n<match.size(); ++n) {
				  // the mapped index of the beginning of the edge
				StStMap::const_iterator fromIdx = m2i.find(match.at(n));
				assert( fromIdx != m2i.end() /* the from-node index has to be available */);
				 // iterate through all adjacent nodes
				for( sgm::Graph_Interface::OutEdge_iterator curEdge = molecules.getOutEdgesBegin(match.at(n)),
						curEdgeEnd = molecules.getOutEdgesEnd(match.at(n));
					curEdge != curEdgeEnd; ++curEdge )
				{
					  // the mapped index of the end of the edge
					StStMap::const_iterator toIdx = m2i.find(curEdge->getToIndex());
					  // check if edge end is among the present nodes 
					  // and if index is higher to avoid duplicated insertions
					  // and if edge was not already present in leftSide
					if (	toIdx != m2i.end() 
							&& curEdge->getToIndex() > curEdge->getFromIndex()
							&& !(boost::edge(fromIdx->second, toIdx->second, its).second)
						) 
					{
						  // add new edge
						CoreGraph::edge_descriptor newEdge
						= boost::add_edge(	boost::vertex( fromIdx->second, its),
											boost::vertex( toIdx->second,   its),
											its ).first;
						  // set context
						itsEdgeContext[newEdge] = ChemRule::RULE_CONTEXT;
						  // set label
						itsEdgeLabel[newEdge] = curEdge->getEdgeLabel();
					}
				}
			}

		}

		  // get first index of each component (direct access)
		const IdxSet & firstOfEach = leftSide.getFirstOfEachComponent();
		const sgm::Graph_Interface::CompLabel & compLabel
				= leftSide.getComponentLabeling();

		  // add the matched graph of each component to its
		for (	IdxSet::const_iterator it = firstOfEach.begin();
				it!=firstOfEach.end(); ++it ) 
		{
			  // if each component should have an own molecule copy, 
			  // add node indices of this component to the mapping
			if (addEachComponent) {
				const size_t curLabel = compLabel.at(*it);
				m2i.clear();
				size_t nodesInComponent = 0;
				  // add all index matches of THIS component
				for (size_t i=0; i<match.size(); ++i ) {
					if (compLabel.at(i) == curLabel) {
						m2i.insert(M2Ival(match.at(i),rule.getLeftSide().nodes.at(i)));
						nodesInComponent++;
					}
				}
				  // check if each node is mapped onto a unique one in target
				if ( m2i.size() != nodesInComponent ) {
					/*************************************************************
					 * some nodes of the leftside rule are mapped onto the same
					 * node in the target graph and no separate copies have to
					 * be created
					 * --> aborting rule application for this match
					 ************************************************************/
					return false;
				}

				// add missing edges from molecule between nodes from the
				// leftside rule pattern, and each only once
				for (StStMap::const_iterator fromIdx=m2i.begin(); fromIdx!=m2i.end(); ++fromIdx) {
					 // iterate through all adjacent nodes
					for( sgm::Graph_Interface::OutEdge_iterator curEdge = molecules.getOutEdgesBegin(fromIdx->first);
						curEdge != molecules.getOutEdgesBegin(fromIdx->first);
						++curEdge )
					{
						  // the mapped index of the end of the edge
						StStMap::const_iterator toIdx = m2i.find(curEdge->getToIndex());
						  // check if edge end is among the present nodes 
						  // and if index is higher to avoid duplicated insertions
						  // and if edge was not already present in rule (==ITS)
						if (	toIdx != m2i.end() 
								&& curEdge->getToIndex() > curEdge->getFromIndex()
								&& !(boost::edge(fromIdx->second, toIdx->second, its).second)
							) 
						{
							  // add new edge
							CoreGraph::edge_descriptor newEdge
							= boost::add_edge(	boost::vertex( fromIdx->second, its),
												boost::vertex( toIdx->second,   its),
												its ).first;
							  // set context
							itsEdgeContext[newEdge] = ChemRule::RULE_CONTEXT;
							  // set label
							itsEdgeLabel[newEdge] = curEdge->getEdgeLabel();
						}
					}
				}

			} // if addEachComponent
			
			  // will contain all indices of molecules that are part of the
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

				  // check if it is a new node
				const bool isNewNode = m2i.find(curNodeIdx) == m2i.end();

				  // add if not already part of the ITS
				if ( isNewNode ) {
					  // add new vertex
					NodeDescr newNode = boost::add_vertex(its);
					  // store mapping
					m2i.insert(M2Ival( curNodeIdx, itsNodeIndex[newNode]));
					  // set context
					itsNodeContext[newNode] = ChemRule::RULE_CONTEXT;
					  // set label
					itsNodeLabel[newNode] = molecules.getNodeLabel(curNodeIdx);
				}
				
				 // iterate through all adjacent nodes
				for( sgm::Graph_Interface::OutEdge_iterator curEdge = molecules.getOutEdgesBegin(curNodeIdx),
						curEdgeEnd = molecules.getOutEdgesEnd(curNodeIdx);
					curEdge != curEdgeEnd; ++curEdge )
				{
					
					  // if neighbored node was not already handled add to
					  // processing queue
					if (handledNodes.find(curEdge->getToIndex())==handledNodes.end()) {
						nodesToHandle.insert(curEdge->getToIndex());
					}
					  // add edge if at least one new node involved and both already added
					if (	m2i.find(curEdge->getToIndex()) != m2i.end()
							&& isNewNode ) 
					{

						  // add new edge
						EdgeDescr newEdge
						= boost::add_edge(	boost::vertex( m2i.find(curEdge->getFromIndex())->second, its),
											boost::vertex( m2i.find(curEdge->getToIndex())->second,   its),
											its).first;
						  // set context
						itsEdgeContext[newEdge] = ChemRule::RULE_CONTEXT;
						  // set label
						itsEdgeLabel[newEdge] = curEdge->getEdgeLabel();
					}
				}
				
			} // while not all node handled
			
		} // for each component
	
		return true;
	}

//##############################################################################

 } // namespace chem
} // namespace ggl

