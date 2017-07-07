#ifndef GGL_MR_APPLYRULE_HH_
#define GGL_MR_APPLYRULE_HH_

#include "sgm/Match_Reporter.hh"
#include "ggl/Graph_Storage.hh"

#include "ggl/Graph.hh"
#include "ggl/RuleGraph.hh"



#include <climits>

namespace ggl {


	  /*! @brief Graph Grammar Rule application for each reported match
	   *
	   *  A sgm::Match_Reporter implementation that applies a graph grammar ggl::Rule
	   *  to a graph and generates the resulting graph. The new graph is added
	   *  to a specified container for further handling like storage or output.
	   *  The pattern graph has to be an instance of ggl::LeftSidePattern.
	   *  
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */ 
	class MR_ApplyRule : public sgm::Match_Reporter {
		
	protected:
		
		  //! the storage to that "reportHit(..)" adds the resulting graphs 
		  //! after the application of a rule to a matched target graph.
		Graph_Storage & storage;
		
		  //! we apply undirected rules or not
		const bool undirectedRule;
		
		
		  /*! if set to true, than all components of the Rules LeftSidePattern 
		   * are matched to an own target graph copy if two components map to 
		   * the same target graph.
		   */
		const bool addEachComponent;
		
	public:
		
		  /*! Construction of a MR_ApplyRule object that adds the resulting
		   *  graphs to the given Graph_Storage object.
		   *  @param storage the Graph_Storage object to add the results to
		   *  @param addEachComponent if set to true, than all components of the
		   *         Rules LeftSidePattern are matched to an own target graph
		   *         copy if two components map to the same target graph.
		   *         NOTE: only important for rules with multi-component 
		   *         LeftSidePattern!
		   */
		MR_ApplyRule(	Graph_Storage & storage
						, const bool addEachComponent = false );
		
		virtual
		~MR_ApplyRule();
		
		
		  //! Applies the Rule represented by the pattern onto the matched 
		  //! subgraph of target and adds the resulting graph to the internal
		  //! Graph_Storage object.
		  //! NOTE: It is assumed that pattern is an instance of 
		  //! ggl::LeftSidePattern!
		  //! @param pattern the pattern graph that was searched for
		  //! @param target the graph the pattern was found within
		  //! @param match contains the indices of the matched pattern nodes in
		  //! the target graph. match[i] corresponds to the mapping of the i-th
		  //! vertex in the pattern graph.
		virtual
		void
		reportHit (	const sgm::Pattern_Interface & pattern,
					const sgm::Graph_Interface & target,
					const sgm::Match & match );


	protected:

		//! Produces the new node label that results from a label rewrite.
		//! This function just returns newNodeLabel but it serves as proxy for
		//! alternative behaviour of node label rewrites
		//! @param oldNodeLabel the old node label (matched by the rule)
		//! @param newNodeLabel the new node label (to be set by the rule)
		//! @return newNodeLabel
		virtual
		std::string
		getAlteredNodeLabel( const std::string & oldNodeLabel, const std::string & newNodeLabel ) const;


	};


////////////////////////////////////////////////////////////////////////////////

} // namespace ggl

#include "ggl/MR_ApplyRule.icc"

#endif /*GGL_MR_APPLYRULE_HH_*/
