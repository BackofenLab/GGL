
#include "ggl/chem/GS_MolCheck.hh"
#include "ggl/chem/AP_NSPDK.hh"
#include "ggl/chem/MR_Reactions.hh"

namespace ggl {
 namespace chem {
	
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



	MR_Reactions::
	MR_Reactions(	const Smiles2GraphMap& smiles2graph_
					, Smiles2GraphMap& newSmiles2graph_
					, Reaction_Container& reactions_
					, const bool addEachComponent_
					, const ReactionRateCalculation * rateCalculation_
					, const AromaticityPerception & aromaticity
					, const bool keepAtomClass
					)
	 :	SuperClass(	newSmiles2graph_
						, smiles2graph_ )
		, graph2smiles()
		, reactions( reactions_ )
		, addEachComponent(addEachComponent_)
		, molchecker( *this, aromaticity )
		, Ruler( keepAtomClass
				? new MR_ApplyRule_NodeLabelPrefix(molchecker, ":", addEachComponent)
				: new MR_ApplyRule( molchecker, addEachComponent ) )
		, rateCalculation(rateCalculation_)
	{
		typedef Graph2SmilesMap::value_type V_T;
		Smiles2GraphMap::const_iterator it;
		  // fill inverse lookup for graph-to-SMILES mapping
		for ( it=smiles2mol->begin(); it != smiles2mol->end(); ++it )
		{
			  // insert inverted iterator for reverse lookup
			graph2smiles.insert( V_T(it->second, it->first) );
		}

		  // fill and check graph-to-SMILES mapping
		bool wasUnknown = true;
		for ( it=smiles2molCheckOnly1->begin();
				it != smiles2molCheckOnly1->end(); ++it )
		{
			  // insert inverted iterator for reverse lookup
			wasUnknown = graph2smiles.insert( V_T(it->second, it->first)).second;
			assert( wasUnknown /* multiple occurrences of SMILES among containers */);
		}
	}


////////////////////////////////////////////////////////////////////////////////

	MR_Reactions::
	MR_Reactions(	const Smiles2GraphMap& smiles2graph_
					, const Smiles2GraphMap& smiles2graph2_
					, Smiles2GraphMap& newSmiles2graph_
					, Reaction_Container& reactions_
					, const bool addEachComponent_
					, const ReactionRateCalculation * rateCalculation_
					, const AromaticityPerception & aromaticity
					, const bool keepAtomClass
					)
	 :	SuperClass(	newSmiles2graph_
						, smiles2graph_
						, smiles2graph2_ )
		, graph2smiles()
		, reactions( reactions_ )
		, addEachComponent(addEachComponent_)
		, molchecker( *this, aromaticity )
		, Ruler( keepAtomClass
				? new MR_ApplyRule_NodeLabelPrefix(molchecker, ":", addEachComponent)
				: new MR_ApplyRule( molchecker, addEachComponent ) )
		, rateCalculation(rateCalculation_)
	{
		typedef Graph2SmilesMap::value_type V_T;
		Smiles2GraphMap::const_iterator it;
		  // fill inverse lookup for graph-to-SMILES mapping
		for ( it=smiles2mol->begin();
				it != smiles2mol->end(); ++it )
		{
			  // insert inverted iterator for reverse lookup
			graph2smiles.insert( V_T(it->second, it->first) );
		}

		bool wasUnknown = false;
		  // fill and check lookup for graph-to-SMILES mapping
		for (it=smiles2molCheckOnly1->begin();
				it != smiles2molCheckOnly1->end(); ++it )
		{
			  // insert inverted iterator for reverse lookup
			wasUnknown = graph2smiles.insert( V_T(it->second, it->first)).second;
			assert( wasUnknown /* multiple occurrences of SMILES among containers */);
		}
		  // fill and check graph-to-SMILES mapping
		for (it=smiles2molCheckOnly2->begin();
				it != smiles2molCheckOnly2->end(); ++it )
		{
			  // insert inverted iterator for reverse lookup
			wasUnknown = graph2smiles.insert( V_T(it->second, it->first)).second;
			assert( wasUnknown /* multiple occurrences of SMILES among containers */);
		}
	}
		

////////////////////////////////////////////////////////////////////////////////
		


	void
	MR_Reactions::
	reportHit (	const sgm::Pattern_Interface & pattern,
				const sgm::Graph_Interface & target,
				const sgm::Match & match )
	{
		
		 // clear content of current reaction
		curReaction.clear();

		  // check for cast success
		assert( dynamic_cast<const LeftSidePattern*>(&pattern) != NULL /* pattern is no LeftSidePattern object ! */);
		  // cast graph to GGL rule
		const LeftSidePattern* leftSide
			= static_cast<const LeftSidePattern*>(&pattern);
		
		/////////////  STORE CHEMRULE ID  ////////////////////////////////////

		  // access to the rule ID that was matched
		curReaction.rule_id = leftSide->getRule().getID();
		
		/////////////  STORE METABOLITES /////////////////////////////////
		
		// the target graph this MR class is able to process
		typedef Molecule_Graph_V TargetGraph;
		// cast graph to processable class
		const TargetGraph* graph
			= dynamic_cast<const TargetGraph*>(&target);

		assert( graph != NULL /* target is no sgm::Graph_boostV_p instance */ );
		
		  // get first index of each component (direct access)
		typedef LeftSidePattern::IndexSet IdxSet;
		const IdxSet & firstOfEach = leftSide->getFirstOfEachComponent();
		
		  // create metabolite SMILES list
		std::set<size_t> matchedGraphIDs;
		  // add SMILES of matched graph of each component to metabolites
		for (	IdxSet::const_iterator it = firstOfEach.begin();
				it!=firstOfEach.end(); ++it )
		{
			  // if not all components are matched onto an own graph copy
			  // check if this molecule was already reported or not
			if (	!addEachComponent
					&& !(matchedGraphIDs.insert(graph->getLocalIndex(match.at( *it )).first).second)
				)
			{
				continue;
			}
			  // get local index for matched position of the first
			  // index in the pattern of the current component
			TargetGraph::LocalIndex idx
				= graph->getLocalIndex(match.at( *it ));
			assert( graph2smiles.find(graph->getGraphs().at(idx.first)) != graph2smiles.end() /* SMILES of mapped molecule is not available in SMILES to graph mapping */);
			  // lookup the corresponding SMILES for the given target graph
			curReaction.metabolites.insert(
					graph2smiles.find(graph->getGraphs().at(idx.first))->second );
		}
		
		/////////////  APPLY CHEMRULE AND STORE PRODUCTS  ////////////////////

		  // forward the hit to the rule applier to collect the resulting
		  // molecules via this->add() function
		Ruler->reportHit( pattern, target, match );
		
		  // check if a molecule was produced, if not abort
		if (curReaction.products.size() == 0) {
			return;
		}
		
		/////////////  ADD REACTION TO STORAGE  //////////////////////////


		  // if a rate calculator is set and the reaction was not known
		  // so far --> calculate reaction rate
		if ( rateCalculation != NULL ) {
			  // try to locate among the known reactions
			Reaction_Container::iterator reactionIt
					= reactions.find( curReaction );
			  // if not already present in container
			if ( reactionIt == reactions.end() ) {
				  // check if transition state is needed
				if (rateCalculation->needTransitionState()) {
					// calculate full imaginary transition state
					ReactionTransitionState its( *leftSide, target, match, addEachComponent );
					  // add SMILES representation to reaction
					curReaction.transState = its.getSMILES();
				}
				 // --> calculate and update reaction rate
				curReaction.rate = rateCalculation->getRate( curReaction );
				 // --> add to container
				reactions.insert( curReaction );
			}
		} else {
			  // --> add to container if not present
			reactions.insert( curReaction );
		}
	}
		

////////////////////////////////////////////////////////////////////////////////
		


	bool
	MR_Reactions::
	insert2map(const std::string & SMILES, const Molecule& graph )
	{
		const bool newSMILES = SuperClass::insert2map( SMILES, graph );
		
		  // if current SMILES was added to storage (unknown so far)
		  // --> add to inverse lookup
		if ( newSMILES ) {
			  // get new entry
			Smiles2GraphMap::iterator newEntry = smiles2mol->find(SMILES);
			  // store inverse information for lookup
#ifndef NDEBUG
			bool wasUnknown =
#endif
					graph2smiles.insert(
						Graph2SmilesMap::value_type( newEntry->second
													, newEntry->first )
#ifndef NDEBUG
					).second;
			assert( wasUnknown /* multiple occurrences of SMILES among containers */);
#else
					);
#endif
		}
		  // store product in current reaction to fill
		curReaction.products.insert(SMILES);
		return newSMILES;
	}


////////////////////////////////////////////////////////////////////////////////

	
	
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 } // namespace chem
} // namespace ggl
