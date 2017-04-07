/*
 * GraphScaffold.cc
 *
 *  Created on: 05.04.2013
 *      Author: mmann
 */

#include "GraphScaffold.hh"

namespace sgm {

	//////////////////////////////////////////////////////////////////////////

	GraphScaffold::GraphScaffold()
	 : RP_Hanser96()
	{
	}

	//////////////////////////////////////////////////////////////////////////

	GraphScaffold::~GraphScaffold() {
	}

	//////////////////////////////////////////////////////////////////////////

	GraphScaffold::
	ScaffoldAnnotation
	GraphScaffold::
	getScaffoldAnnotation( const Graph_Interface & graph
							, const size_t maxRingSize )
	{
		return getScaffoldAnnotation( graph, NULL, maxRingSize );
	}

	//////////////////////////////////////////////////////////////////////////

	GraphScaffold::
	ScaffoldAnnotation
	GraphScaffold::
	getScaffoldAnnotation( const Graph_Interface & graph
							, RingReporter & nextReporter
							, const size_t maxRingSize )
	{
		return getScaffoldAnnotation( graph, &nextReporter, maxRingSize );
	}

	//////////////////////////////////////////////////////////////////////////

	GraphScaffold::
	ScaffoldAnnotation
	GraphScaffold::
	getScaffoldAnnotation( const Graph_Interface & graph
							, RingReporter * nextReporter
							, const size_t maxRingSize )
	{
		  // initialize annotation
		ScaffoldAnnotation annotation(graph.getNodeNumber(),GST_UNKNOWN);

		  // create ring annotation reporter
		RR_Annotation reporter( annotation, nextReporter );

		  // create pGraph of the given graph
		initializeP_Graph( graph, reporter );

		bool danglingsLeft = true;

		  // iterative collapse of pGraph
		while( ! toRemove.empty() ) {

			  // sort remaining indices to be removed by decreasing degree
			std::sort( toRemove.begin(), toRemove.end(), degree_sort(pGraphDegree));

			danglingsLeft = danglingsLeft && pGraphDegree.at(*(toRemove.rbegin())) <= 1;

			  // check for not yet annotated dangling ends
			if ( danglingsLeft
					&& annotation.at(*(toRemove.rbegin())) == GST_UNKNOWN)
			{
				  // set annotation accordingly
				annotation[*(toRemove.rbegin())] = GST_DANGLING;
			}

			  // remove all edges from the given vertex and report rings if any
			remove_vertex( *(toRemove.rbegin()), graph, reporter, maxRingSize );

			  // remove the handled entry (last entry)
			toRemove.resize( toRemove.size()-1 );
		}

		  // set all remaining nodes with unknown type to linker
		for (ScaffoldAnnotation::iterator a = annotation.begin(); a != annotation.end(); ++a) {
			if (*a == GST_UNKNOWN) {
				*a = GST_LINKER;
			}
		}

		  // return final annotation
		return annotation;
	}

	//////////////////////////////////////////////////////////////////////////

	void
	GraphScaffold::
	RR_Annotation::
	reportRing( const Graph_Interface& graph, const RingList & ringList )
	{
		  // store ring annotation
		for (RingList::const_iterator r = ringList.begin(); r != ringList.end(); ++r ){
			annotation[*r] = GST_RING;
		}

		  // forward ring report call
		if (nextReporter != NULL) {
			nextReporter->reportRing( graph, ringList );
		}
	}

	//////////////////////////////////////////////////////////////////////////


} // namespace sgm
