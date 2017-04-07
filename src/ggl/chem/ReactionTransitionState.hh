#ifndef GGL_CHEM_REACTIONTRANSITIONSTATE_HH_
#define GGL_CHEM_REACTIONTRANSITIONSTATE_HH_

#include "ggl/chem/ChemRule.hh"
#include "ggl/chem/ChemRuleGraph.hh"

#include <sgm/Graph_Interface.hh>
#include <sgm/Match_Reporter.hh>


namespace ggl {
 namespace chem {


	/*! @brief Transition state of a reaction
	 *
	 * This class encodes an imaginary transition state (ITS) of a chemical
	 * reaction. 
	 * 
	 * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	class ReactionTransitionState : public ChemRule::TransitionState
	{
	public:
		
		  //! access to the superclass TransitionState type
		typedef ChemRule::TransitionState SuperClass;
		
		
	protected:
		
		  //! the graph encoding of the imaginary transition state
		using SuperClass::tState;
		
		  //! SMILES representation of the imaginary transition state
		using SuperClass::tStateSMILES;
		
	public:
		
		 /*! Constructs an imaginary transition state based on the educts, 
		  * products, and chemical GGL rule applied.
		  * 
		  * @param leftSide the left side pattern of the chemical rule applied
		  * @param molecules the molecules where leftSide was mapped on
		  * @param match the mapping of nodes of leftSide on the nodes of the
		  *         molecules object
		  * @param addEachComponent if set to true, than all components of the
		  *        Rules LeftSidePattern are matched to an own molecule graph
		  *        copy if two components map to the same molecule graph.
		  */
		ReactionTransitionState(
					const LeftSidePattern & leftSide
					, const sgm::Graph_Interface & molecules
					, const sgm::Match & match
					, const bool addEachComponent = false
				);
		
		virtual ~ReactionTransitionState();
		
	protected:
		
		 /*! Merges the given molecules into a rule graph such that the 
		  * left side pattern of the rule maps the matched subgraph of the 
		  * molecules object. This is the first step to create the imaginary
		  * transition state.
		  * 
		  * @param its the intermediate transition state object to extend
		  *        (NOTE : the object will be overwritten during the merge!)
		  * @param leftSide the left side pattern of the chemical rule applied
		  * @param molecules the molecule information to add to 
		  * @param match the mapping of nodes of leftSide on the nodes of the
		  *         molecules object
		  * @param addEachComponent if set to true, than all components of the
		  *        Rules LeftSidePattern are matched to an own molecule graph
		  *        copy if two components map to the same molecule graph.
		  */
		static
		bool
		mergeGraphs(	ChemRule::CoreGraph & its
						, const LeftSidePattern & leftSide
						, const sgm::Graph_Interface & molecules
						, const sgm::Match & match
						, const bool addEachComponent = false );
	};
	
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 } // namespace chem
} // namespace ggl

// implementation
#include "ggl/chem/ReactionTransitionState.icc"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	
#endif /*REACTIONTRANSITIONSTATE_HH_*/

