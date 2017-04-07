/*
 * VF2_MatchingHandler.cc
 *
 *  Created on: 19.11.2011
 *      Author: mmann
 */

#include "VF2_MatchingHandler.hh"

#include <limits>
#include <sstream>
#include <map>

namespace sgm {

////////////////////////////////////////////////////////////////////////////////

	VF2_MatchingHandler
	::VF2_MatchingHandler()
	{
	}

////////////////////////////////////////////////////////////////////////////////

	VF2_MatchingHandler
	::~VF2_MatchingHandler()
	{
	}

////////////////////////////////////////////////////////////////////////////////


	bool
	VF2_MatchingHandler
	::fillLoader(	const Graph_Interface& graph,
					vf2::ARGEdit & loader,
					PatternMap & graphLabel,
					NodeDataVec & nodeData,
					EdgeLabelVec & edgeLabels,
					bool isPattern )
	{
		  // container for edge label handling
		typedef
#if HAVE_UNORDERED_MAP > 0
			std::unordered_map<Graph_Interface::IndexType, EdgeLabel>
#elif HAVE_TR1_UNORDERED_MAP > 0
			std::tr1::unordered_map<Graph_Interface::IndexType, EdgeLabel>
#elif HAVE_GNU_HASH_MAP > 0
			__gnu_cxx::hash_map<Graph_Interface::IndexType, EdgeLabel>
#else
			std::map<Graph_Interface::IndexType, EdgeLabel>
#endif
				EdgeMap;

		  // mapping of EdgeLabels on indices within edgeLabels to enable a
		  // reuse of edgelabel objects when filling the loader
		typedef std::map< EdgeLabel, size_t > EdgeLabel2idxMap;

		bool noAdditionalLabels = true;

		const size_t nodeNumber = graph.getNodeNumber();

		  // resize node data container
		nodeData.resize(nodeNumber);
		 // add all nodes
		for (size_t i=0; i<nodeNumber; ++i) {
			std::string label = graph.getNodeLabel(i);
			  // check if label already exists
			PatternMap::iterator pos = graphLabel.find(label);
			  // if not add to map
			if (pos == graphLabel.end()) {
				noAdditionalLabels = false;
				if (isPattern) {
					std::pair<PatternMap::iterator, bool>
					res = graphLabel.insert(PatternMap::value_type(label, graphLabel.size()+1 ));
					pos = res.first;
				}
			}
			  // store node data object
			nodeData[i] = NodeData(i,pos == graphLabel.end() ? 0 : pos->second, 0, 0 );
			  // add labeled node to loader
			loader.InsertNode( & nodeData[i] );
		}

		  // store known edge labels
		EdgeLabel2idxMap edgeLabel2idx;
		for (size_t i=0; i<edgeLabels.size(); ++i) {
			edgeLabel2idx[ *edgeLabels.at(i) ] = i;
		}

		 // add all edges
		for (size_t i=0; i<nodeNumber; ++i) {
			EdgeMap edges;
			EdgeMap::iterator edgeIt;
			size_t outDegree = 0;
			size_t selfloops = 0;
			  // collect parallel edges
			for (	Graph_Interface::OutEdge_iterator e = graph.getOutEdgesBegin(i),
					eEnd = graph.getOutEdgesEnd(i); e != eEnd; ++e )
			{
				if (e->getFromIndex() != e->getToIndex()) {
					++outDegree;
				} else {
					++selfloops;
				}
				  // get access to edge label instance or create if non-existent
				if (edges.find(e->getToIndex()) == edges.end()) {
					edges[e->getToIndex()] = EdgeLabel();
				}
				EdgeLabel &curLabels = edges.find(e->getToIndex())->second;

				  // check if label already exists
				PatternMap::iterator pos = graphLabel.find(e->getEdgeLabel());
				  // if not add to map
				if (pos == graphLabel.end()) {
					noAdditionalLabels = false;
					if (isPattern) {
						pos = graphLabel.insert(PatternMap::value_type(e->getEdgeLabel(), graphLabel.size()+1 )).first;
					}
				}
				  // add current out-label to labels
				curLabels.insert(pos == graphLabel.end() ? 0 : pos->second);
			}
			  // set node out degree information
			nodeData[i].outDegree = outDegree;
			  // set self loop information
			nodeData[i].selfloops = selfloops;
			  // add collected edges to loader
			EdgeLabel2idxMap::const_iterator edgeLabelIdx;
			for (EdgeMap::iterator e=edges.begin(); e!=edges.end(); ++e) {
				  // get edge label
				size_t edgeLabelIdx = 0;

				  // check if edge label known
				if (edgeLabel2idx.find( e->second ) == edgeLabel2idx.end()) {
					  // not known yet -> store
					edgeLabelIdx = edgeLabels.size();
					edgeLabels.push_back( new EdgeLabel(e->second) );
					edgeLabel2idx[ e->second ] = edgeLabelIdx;
				} else {
					edgeLabelIdx = edgeLabel2idx.find( e->second )->second;
				}
				  // add edge to loader
				loader.InsertEdge(	i,
									e->first,
									edgeLabels.at(edgeLabelIdx) );
			}
		}

		  /////////////////////  STORE CONSTRAINTS  //////////////////////////

		if (isPattern) {
			  // get access to pattern functionality
			const Pattern_Interface * pattern = dynamic_cast<const Pattern_Interface*>(&graph);
			  // check if this is really a pattern
			if ( pattern != NULL ) {
				  // get constraints
				const Pattern_Interface::ConstraintVec & constraints = pattern->getConstraints();
				const MC_Node * nodeConstraint;
				for (size_t i=constraints.size(); i>0; --i) {
					  // check if current constraint is a node constraint
					nodeConstraint = dynamic_cast<const MC_Node*>(constraints.at(i-1));
					if (nodeConstraint != NULL) {
						  // check if we have already constraints for this node
						assert(nodeConstraint->constrainedNodeID<nodeData.size());
						if (nodeData.at(nodeConstraint->constrainedNodeID).nodeConstraints == NULL) {
							  // create container
							nodeData[nodeConstraint->constrainedNodeID].nodeConstraints
								= new NodeData::NodeConstraints();
						}
						  // store constraint in according constraint container
						nodeData[nodeConstraint->constrainedNodeID].nodeConstraints->push_back(nodeConstraint);
						  // proceed to next constraint
						continue;
					}
					  // TODO add edge constraint handling
				}
			}

		}


		  // whether or not additional labels have been found
		return isPattern || noAdditionalLabels;
	}

////////////////////////////////////////////////////////////////////////////////


	bool
	VF2_MatchingHandler
	::isApplicable (	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						const Match & match )
	{

		 // access to the additional matching constraints
		const Pattern_Interface::ConstraintVec & mC = pattern.getConstraints();

		  // fast check if any constraint present or not
		if ( mC.empty() ) {
			return true;
		}

		bool allConstraintsMatched = true;

		  // check all additional constraints
		for (sgm::Pattern_Interface::ConstraintVec::const_iterator rC = mC.begin();
				allConstraintsMatched && rC != mC.end(); ++rC )
		{
			  // check if current constraint is violated
			allConstraintsMatched = (*rC)->isValidMatch( pattern, target, match );
		}

		return allConstraintsMatched;
	}

////////////////////////////////////////////////////////////////////////////


	bool
	VF2_MatchingHandler
	::areDegreeCompatible(	const VF2_MatchingHandler::NodeDataVec & patternNodes
						, const VF2_MatchingHandler::NodeDataVec & targetNodes )
	{
		typedef
#if HAVE_UNORDERED_MAP > 0
			std::unordered_map<int, size_t>
#elif HAVE_TR1_UNORDERED_MAP > 0
			std::tr1::unordered_map<int, size_t>
#elif HAVE_GNU_HASH_MAP > 0
			__gnu_cxx::hash_map<int, size_t>
#else
			std::map<int, size_t>
#endif
			DegreeDistribution;

		DegreeDistribution patternDegreeCount;

		int maxDegree = -1;
		int minDegree = INT_MAX;

		size_t nodes2match = patternNodes.size();

		  // gather pattern node degree information
		for (NodeDataVec::const_iterator p=patternNodes.begin(); p!=patternNodes.end(); ++p) {
			if (patternDegreeCount.find((int)(p)->outDegree)==patternDegreeCount.end())
				patternDegreeCount[ (int)(p)->outDegree ] = 1;
			else
				patternDegreeCount[ (int)(p)->outDegree ] += 1;
			  // check for max degree
			if ( maxDegree < (int)(p)->outDegree )
				maxDegree = (int)(p)->outDegree;
			  // check for min degree
			if ( minDegree > (int)(p)->outDegree )
				minDegree = (int)(p)->outDegree;
		}

		  // check target if it shows at least the pattern degree distribution
		for (NodeDataVec::const_iterator t=targetNodes.begin();
				nodes2match > 0 && maxDegree >= minDegree && t!=targetNodes.end(); ++t)
		{
			  // node with more than enough out edges
			if ((int)(t)->outDegree >= maxDegree) {
				patternDegreeCount[maxDegree] -= 1;
				  // decrease matching counter
				--nodes2match;
				  // check if all nodes of maximal degree are matched
				if (patternDegreeCount[maxDegree] == 0) {
					--maxDegree;
					  // find next smaller populated degree
					while( maxDegree >= minDegree
							&& (patternDegreeCount.find(maxDegree) == patternDegreeCount.end()
									|| patternDegreeCount[maxDegree] == 0) )
					{ --maxDegree; }
				}
			}
			  // find next equal or smaller pattern degree to decrease counter
			else if ((int)(t)->outDegree >= minDegree){
				  // find degree counter to decrease
				int curDegree = (int)(t)->outDegree;
				while(curDegree >= minDegree && (patternDegreeCount.find(curDegree) == patternDegreeCount.end()
									|| patternDegreeCount[curDegree] == 0))
				{ --curDegree; }
				  // decrease degree counter
				patternDegreeCount[curDegree] -= 1;
				  // decrease matching counter
				--nodes2match;
				  // check if minimal degree has to be updated
				if (curDegree == minDegree && patternDegreeCount[minDegree] == 0) {
					++minDegree;
					  // find next populated degree level
					while(minDegree <= maxDegree && (patternDegreeCount.find(minDegree) == patternDegreeCount.end()
							|| patternDegreeCount[minDegree] == 0))
					{ ++minDegree; }
				}
			}
			  // else ignore node since uninteresting degree
		}

		return nodes2match == 0;

	}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	VF2_MatchingHandler
	::VF2_match_handler
	::VF2_match_handler( Match_Reporter & matchReporter
						, const Pattern_Interface& pattern
						, const Graph_Interface& target
						, const size_t maxHitsToFind
						)
	 :	matchReporter(matchReporter)
		, pattern(pattern)
		, target(target)
		, maxHitsToFind(maxHitsToFind)
		, numOfMatches(0)
	{
	}

////////////////////////////////////////////////////////////////////////////////

	bool
	VF2_MatchingHandler
	::VF2_match_handler
	::report_matches (int n, vf2::node_id *q, vf2::node_id *t, void *x)
	{

		  // cast user data into data object
		VF2_match_handler * data = static_cast<VF2_match_handler*>(x);

		  // fill reporter data structure
		Match res(n,0);
		for(int i=0; i<n; ++i) {
			res[(size_t)q[i]] = (Graph_Interface::IndexType)t[i];
		}
		  // check if additional constraints are fulfilled
		if (VF2_MatchingHandler::isApplicable( data->pattern, data->target, res)) {
			  // report match
			data->matchReporter.reportHit( data->pattern, data->target, res);
			  // count match
			data->numOfMatches += 1;
		}
		  // check if we stop the computation since enough matches have been found
		if (data->numOfMatches >= data->maxHitsToFind) {
			return true;
		}

		  // false == dont stop computation
		return false;
	}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


	VF2_MatchingHandler
	::NodeData
	::NodeData( )
	 :	id(std::numeric_limits<size_t>::max())
		, label(0)
		, outDegree(0)
		, selfloops(0)
		, nodeConstraints(NULL)
	{}



////////////////////////////////////////////////////////////////////////////////


	VF2_MatchingHandler
	::NodeData
	::NodeData( const size_t id
				, VF2_MatchingHandler::Label label
				, const size_t outDegree
				, const size_t selfloops )
	 :	id(id)
		, label(label)
		, outDegree(outDegree)
		, selfloops(selfloops)
		, nodeConstraints(NULL)
	{}


////////////////////////////////////////////////////////////////////////////////


	VF2_MatchingHandler
	::NodeData
	::NodeData( const VF2_MatchingHandler::NodeData & toCopy )
	 :	id(toCopy.id)
	 	, label(toCopy.label)
		, outDegree(toCopy.outDegree)
		, selfloops(toCopy.selfloops)
		, nodeConstraints(NULL)
	{
		if (toCopy.nodeConstraints != NULL) {
			  // create a copy
			nodeConstraints = new NodeConstraints(*(toCopy.nodeConstraints));
		}
	}


////////////////////////////////////////////////////////////////////////////////


	VF2_MatchingHandler
	::NodeData
	::~NodeData()
	{
		if (nodeConstraints != NULL) {
			nodeConstraints->clear();
			delete nodeConstraints;
		}
	}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


	VF2_MatchingHandler
	::LabelComparator
	::LabelComparator( const Label pWildcard )
	 :	pWildcard(pWildcard)
	{}


////////////////////////////////////////////////////////////////////////////////


	bool
	VF2_MatchingHandler
	::LabelComparator
	::compatibleLabel(const Label& a, const Label& b) {
		  // check for wildcards
		if (pWildcard != 0 && (a == pWildcard || b == pWildcard)) {
			return true;
		}
		  // if one label is NULL, it was not present in the other graph
		  // i.e. no match
		if (a == 0 || b == 0) {
			return false;
		}
		  // check for equality
		return a == b;
	}


////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}
