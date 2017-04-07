
#include <algorithm>

#include "sgm/RP_Hanser96.hh"

namespace sgm {

//////////////////////////////////////////////////////////////////////////////

RP_Hanser96
::RP_Hanser96()
 :	pGraph()
	, pGraphIndex()
	, pGraphPath()
	, pGraphDegree()
	, toRemove()
{
}

//////////////////////////////////////////////////////////////////////////////


RP_Hanser96
::~RP_Hanser96()
{
}

//////////////////////////////////////////////////////////////////////////////

size_t
RP_Hanser96
::findRings( const Graph_Interface & graph, RingReporter & reporter, const size_t maxRingSize )
{
	  // for counting all reported rings
	size_t ringCount = 0;

	  // create pGraph of the given graph
	ringCount += initializeP_Graph( graph, reporter );

	  // iterative collapse of pGraph
	while( ! toRemove.empty() ) {

		  // sort remaining indices to be removed by decreasing degree
		std::sort( toRemove.begin(), toRemove.end(), degree_sort(pGraphDegree));

		  // remove all edges from the given vertex and report rings if any
		ringCount += remove_vertex( *(toRemove.rbegin()), graph, reporter, maxRingSize );

		  // remove the handled entry (last entry)
		toRemove.resize( toRemove.size()-1 );
	}

	  // final number of rings found
	return ringCount;
}


//////////////////////////////////////////////////////////////////////////////

size_t
RP_Hanser96::
findRingBonds( const Graph_Interface & graph
			, BondSet & ringBonds
			, const size_t maxRingSize )
{
	  // for counting all reported rings
	ringBonds.clear();

	// TODO : there are better implementation possibilities ...

	  // create pGraph of the given graph and store single bond loops
	LoopBondReporter reporter(ringBonds);
	findRings( graph, reporter, maxRingSize );

	  // get final number of ring bonds
	return ringBonds.size();
}

//////////////////////////////////////////////////////////////////////////////


size_t
RP_Hanser96
::initializeP_Graph( const Graph_Interface & graph
					, RingReporter& reporter )
{
	size_t ringCount = 0;

	  // rewrite pGraph
	pGraph = P_Graph();
	  // reset access
	pGraphIndex = boost::get( boost::vertex_index_t(), pGraph );
	pGraphPath = boost::get( boost::edge_name_t(), pGraph );

	  // reset data structures
	pGraphDegree.resize(graph.getNodeNumber(),0);
	toRemove.resize(graph.getNodeNumber(),0);

	  // add all vertices
	for (size_t i=0; i<graph.getNodeNumber(); ++i) {
		boost::add_vertex( pGraph );
		pGraphDegree[i] = 0;
		toRemove[i] = i;
	}
	  // add all edges
	for (size_t i=0; i<graph.getNodeNumber(); ++i) {
		for (sgm::Graph_Interface::OutEdge_iterator e = graph.getOutEdgesBegin(i),
				eEnd = graph.getOutEdgesEnd(i); e != eEnd; ++e )
		{
			  // handle direct loops on a single node
			if ( e->getFromIndex() == e->getToIndex()) {
				RingReporter::RingList curRing;
				curRing.push_back( e->getFromIndex() );
				reporter.reportRing( graph, curRing );
				ringCount++;
			} else
			  // handle all other edges but each only once
			  // (therefore order check)
			if ( e->getFromIndex() < e->getToIndex()) {
				  // gather ring information for current edge
				RingInfo curRing;
				curRing.nodes.insert( e->getFromIndex() );
				curRing.nodes.insert( e->getToIndex() );
				curRing.path.push_back( e->getFromIndex() );
				curRing.path.push_back( e->getToIndex() );

				  // insert new edge
				std::pair< P_Graph::edge_descriptor, bool>
				newEdge = boost::add_edge(	boost::vertex( e->getFromIndex(), pGraph )
											, boost::vertex( e->getToIndex(), pGraph )
											, pGraph );

				  // set edge label = path information
				pGraphPath[ newEdge.first ] = curRing;
				  // update degree information
				pGraphDegree[ e->getFromIndex() ] += 1;
				pGraphDegree[ e->getToIndex() ] += 1;
			}
		}
	}
	  // number of self loops within the input graph
	return ringCount;
}


//////////////////////////////////////////////////////////////////////////////


size_t
RP_Hanser96
::remove_vertex(	const size_t nextToRemove
					, const Graph_Interface& graph
					, RingReporter& reporter
					, const size_t maxRingSize )
{
	size_t ringCount = 0;

	  // access to the neighbored nodes
	P_Graph::out_edge_iterator n1, n1_last, n2, n2_last;

	  // check if path merge is needed
	if ( pGraphDegree.at(nextToRemove) > 1 ) {

		  // get access to all pairs of out edges (no ring among them anymore)
		for( boost::tie( n1, n1_last ) = boost::out_edges( boost::vertex( nextToRemove, pGraph ), pGraph );
				n1 != n1_last; ++n1 )
		{
			boost::tie( n2, n2_last ) = boost::out_edges( boost::vertex( nextToRemove, pGraph ), pGraph );
			while( n2 != n2_last && n2 != n1 ) {
				n2++;
			}
			  // check if end reached
			if ( n2 == n2_last ) {
				continue;
			}
			for( n2++; n2 != n2_last; ++n2 )
			{
				  // add new edge with merged Ring
				const RingInfo& p1 = pGraphPath[ *n1 ];
				const RingInfo& p2 = pGraphPath[ *n2 ];
				size_t curPathFrom = boost::target( *n1, pGraph);
				size_t curPathTo = boost::target( *n2, pGraph);
				  // create joined path information
				RingInfo curPath;
				curPath.nodes.insert(p1.nodes.begin(),p1.nodes.end());
				curPath.nodes.insert(p2.nodes.begin(),p2.nodes.end());
				  // check if intersection of p1 and p2 is only current index or
				  // we have a ring closure
				if (	(p1.nodes.size()+p2.nodes.size()-curPath.nodes.size()) > 2
						|| ((p1.nodes.size()+p2.nodes.size()-curPath.nodes.size()) == 2
								&& curPathFrom != curPathTo) )
				{
					  // the overlap of p1 and p2 is more than the path ends..
					  // --> cannot be a ring anymore
					continue;
				}

				  // create joined path information
				if (*p1.path.begin() == *p2.path.begin()) {
					curPath.path.insert( curPath.path.end(),p1.path.rbegin(), p1.path.rend() );
					curPath.path.insert( curPath.path.end(),++(p2.path.begin()), p2.path.end() );
				} else if ( *p1.path.rbegin() == *p2.path.rbegin() ) {
					curPath.path.insert( curPath.path.end(),p1.path.begin(), p1.path.end() );
					curPath.path.insert( curPath.path.end(),++(p2.path.rbegin()), p2.path.rend() );
				} else if ( *p1.path.rbegin() == *p2.path.begin() ) {
					curPath.path.insert( curPath.path.end(),p1.path.begin(), p1.path.end() );
					curPath.path.insert( curPath.path.end(),++(p2.path.begin()), p2.path.end() );
				} else if ( *p1.path.begin() == *p2.path.rbegin() ) {
					curPath.path.insert( curPath.path.end(),p1.path.rbegin(), p1.path.rend() );
					curPath.path.insert( curPath.path.end(),++(p2.path.rbegin()), p2.path.rend() );
				}

				// proceed only for curRing fragments that do NOT EXCEED maxRingSize
				if (curPath.nodes.size() <= maxRingSize) {

					  // check if ring already formed
					if (curPathFrom == curPathTo) {
						  // report and count
						reporter.reportRing( graph, curPath.path );
						ringCount++;
					} else {
					  // no ring was formed, curPath is really a path

						  // --> create new edge
						std::pair< P_Graph::edge_descriptor, bool>
						newEdge = boost::add_edge(	boost::vertex( curPathFrom, pGraph )
													, boost::vertex( curPathTo, pGraph )
													, pGraph );
						  // set path information
						pGraphPath[ newEdge.first ] = curPath;
						  // update degree information
						pGraphDegree[ boost::vertex( curPathFrom, pGraph ) ] += 1;
						pGraphDegree[ boost::vertex( curPathTo, pGraph ) ] += 1;
					}
				}
			}
		}
	}
	  // remove all adjacent edges of current node
	boost::tie( n1, n1_last ) = boost::out_edges( boost::vertex( nextToRemove, pGraph ), pGraph );
	while ( n1 != n1_last ) {
		  // update degree information
		pGraphDegree[ boost::source( *n1, pGraph ) ] -= 1;
		pGraphDegree[ boost::target( *n1, pGraph ) ] -= 1;
		  // remove edge
		boost::remove_edge( *n1, pGraph );
		  // update iterators
		boost::tie( n1, n1_last ) = boost::out_edges( boost::vertex( nextToRemove, pGraph ), pGraph );
	}

	  // final ring count of this removal iteration
	return ringCount;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

RP_Hanser96::
LoopBondReporter::
LoopBondReporter( BondSet & bondSet )
 : bondSet(bondSet)
{
}

//////////////////////////////////////////////////////////////////////////////
RP_Hanser96::
LoopBondReporter::
~LoopBondReporter()
{
}

//////////////////////////////////////////////////////////////////////////////


void
RP_Hanser96::
LoopBondReporter::
reportRing( const Graph_Interface& graph, const RingList & ringList)
{

	// check if this is a ring at all
	if (ringList.empty()) {
		return;
	}
	// check for loop
	if (ringList.size() == 1) {
		bondSet.insert( BondSet::value_type( *(ringList.begin()), *(ringList.begin()) ) );
		return;
	}
	// at least of length 2;
	RingList::const_iterator source=ringList.begin(), target=ringList.begin();
	for (++target; target != ringList.end(); ++source, ++target) {
		  // add in increasing index order to bond list
		if (*source < *target) {
			bondSet.insert( BondSet::value_type( *source, *target ) );
		} else {
			bondSet.insert( BondSet::value_type( *target, *source ) );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////


} // namespace sgm
