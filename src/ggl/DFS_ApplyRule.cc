
#include "ggl/DFS_ApplyRule.hh"

#include "sgm/MR_SymmBreak.hh"

#include "ggl/MR_ApplyRule.hh"

namespace ggl {

	//////////////////////////////////////////////////////////////////////////


	bool
	DFS_ApplyRule
	::findSolution( const std::vector< Rule > & rules
					, const Graph & startGraph
					, Graph_Storage & solutionStorage_
					, DFS_ApplyRule::DFS_Visitor & visitor_
					, DFS_ApplyRule::SearchTrace * finalTrace_
					, const bool doSymmBreak
					, sgm::SubGraphMatching * sgm_ )
	{
		  // setup rule patterns for matching
		rulePattern.resize(rules.size(), NULL);
		  // setup symmetry breaker container to be filled if needed
		std::vector< sgm::MR_SymmBreak > symmBreaker;
		  // set up Rule applier that reports to this for DFS handling
		MR_ApplyRule mr_applyRule( *this );


		  // initialize reporter with direct rule application without symm break
		reporter.clear();
		reporter.resize( rules.size(), &mr_applyRule );

		  // setup rule pattern and reporting
		for (size_t i=0; i<rules.size(); ++i ) {
			  // create left side pattern
			ggl::LeftSidePattern * pat = new ggl::LeftSidePattern(rules[i]);
			assert(pat != NULL);
			  // store pattern for matching
			rulePattern[i] = pat;

			  // create match reporter
			if (doSymmBreak) {
				  // create and store new symmetry breaker
				symmBreaker.push_back( sgm::MR_SymmBreak( pat->getGraphAutomorphism(), mr_applyRule) );
			}
		}
		  // setup correct symmetry breaking match reporter if needed
		for (size_t i=0; doSymmBreak && i<symmBreaker.size(); ++i ) {
			reporter[i] = & symmBreaker.at(i);
		}

		  // setup DFS variables
		solutionStorage = &solutionStorage_;
		visitor = &visitor_;
		  // create default subgraph matcher if needed
		sgm = sgm_;
		if (sgm == NULL) {
			sgm = new sgm::SGM_vf2();
		}
		finalTrace = finalTrace_;
		currentTrace.clear();


		///////////////////  DFS START /////////////////////////////////////

		  // flag that will hold the search success status
		solutionFound = false;
		try {
			  // start DFS search
			add( startGraph );

		} catch ( DFS_ABORTION & solutionTrigger ) {
			// catch search abortion
		}
		////////////////////  DFS END //////////////////////////////////////


		  // garbage collection
		for (size_t i=0; i<rulePattern.size(); ++i ) {
			delete rulePattern[i];
			rulePattern[i] = NULL;
		}
		rulePattern.clear();
		reporter.clear();
		solutionStorage = NULL;
		visitor = NULL;
		if (sgm != NULL && sgm != sgm_) {
			delete sgm;
		}
		sgm = NULL;
		finalTrace = NULL;
		currentTrace.clear();
		symmBreaker.clear();

		return solutionFound;
	}

	//////////////////////////////////////////////////////////////////////////


	void
	DFS_ApplyRule
	::add( const Graph & graph )
	{
		assert(visitor != NULL);

		  // get DFS status
		DFS_Visitor::Decision lastDecision = visitor->status(graph);

		switch(lastDecision) {

		  // check if search failure and trace back needed
		case DFS_Visitor::FAILURE_STOP :
		case DFS_Visitor::FAILURE_TRACEBACK :
		{
			  // stop search via pseudo exception
			if (lastDecision == DFS_Visitor::FAILURE_STOP) {
				throw DFS_ABORTION();
			}
			return;
		}

		  // check if solution found -> start solution handling
		case DFS_Visitor::SOLUTION_STOP :
		case DFS_Visitor::SOLUTION_TRACEBACK :
		{

			  // store solution
			assert( solutionStorage != NULL );
			solutionStorage->add(graph);
			solutionFound = true;

			  // store trace if needed
			if (finalTrace != NULL) {
				  // clear trace
				for ( size_t i=0; i<finalTrace->size(); ++i ) {
					delete (*finalTrace)[i];
				}
				  // copy trace;
				finalTrace->resize( currentTrace.size()+1, NULL );
				for ( size_t i=0; i<currentTrace.size(); ++i ) {
					(*finalTrace)[i] = new Graph( *(currentTrace.at(i)) );
				}
				  // add final solution
				(*finalTrace)[currentTrace.size()] = new Graph( graph );
			}

			  // stop search via pseudo exception
			if (lastDecision == DFS_Visitor::SOLUTION_STOP) {
				throw DFS_ABORTION();
			}
			return;
		}

		  // check if another DFS iteration is to be done
		case DFS_Visitor::CONTINUE : {

			  // create new target graph
			const Graph_boost target(graph);

			  // add current graph to trace
			currentTrace.push_back( & graph );

			  // apply rules
			assert(sgm != NULL);
			sgm->findMatches( rulePattern, target, reporter, UINT_MAX );

			  // remove current graph from trace
			currentTrace.pop_back( );
			return;
		}

		default :
			assert(false) /*should never happen */;
		}
	}

	//////////////////////////////////////////////////////////////////////////


}  // namespace ggl
