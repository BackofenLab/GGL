#ifndef GGL_MR_APPLYRULE_NODELABELPREFIX_HH_
#define GGL_MR_APPLYRULE_NODELABELPREFIX_HH_


#include "ggl/MR_ApplyRule.hh"


#include <climits>

namespace ggl {


	  /*! @brief Graph Grammar Rule application for each reported match that
	   * applies node label changes only to prefixes
	   *
	   *  A ggl::MR_ApplyRule variant that applies a graph grammar ggl::Rule
	   *  to a graph and generates the resulting graph but applies node label
	   *  changes only to prefixes of the original node labels.
	   *  Thus, suffixes are preserved. The new graph is added
	   *  to a specified container for further handling like storage or output.
	   *  The pattern graph has to be an instance of ggl::LeftSidePattern.
	   *  
	   *  @author Martin Mann (c) 2017 http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */ 
	class MR_ApplyRule_NodeLabelPrefix : public ggl::MR_ApplyRule {
		
	protected:
		
		  /*! the separator that marks the beginning of the node label part to
		   * be preserved by node label rewrites
		   */
		const std::string nodeLabelSeparator;
		
	public:
		
		  /*! Construction of a MR_ApplyRule_NodeLabelPrefix object that adds the resulting
		   *  graphs to the given Graph_Storage object.
		   *  @param storage the Graph_Storage object to add the results to
		   *  @param nodeLabelSeparator the separator that marks the beginning
		   *         of the node label suffix to be preserved from rewrite
		   *  @param addEachComponent if set to true, than all components of the
		   *         Rules LeftSidePattern are matched to an own target graph
		   *         copy if two components map to the same target graph.
		   *         NOTE: only important for rules with multi-component 
		   *         LeftSidePattern!
		   */
		MR_ApplyRule_NodeLabelPrefix(	Graph_Storage & storage
						, const std::string & nodeLabelSeparator
						, const bool addEachComponent = false );
		
		virtual
		~MR_ApplyRule_NodeLabelPrefix();


	protected:

		//! Produces the new node label that preserves the old node label suffix
		//! beginning at a given separator nodeLabelSeparator.
		//! @param oldNodeLabel the old node label (matched by the rule)
		//! @param newNodeLabel the new node label (to be set by the rule)
		//! @return newNodeLabel + suffix_of_oldNodeLabel
		virtual
		std::string
		getAlteredNodeLabel( const std::string & oldNodeLabel, const std::string & newNodeLabel ) const;


	};


////////////////////////////////////////////////////////////////////////////////

} // namespace ggl

#include "ggl/MR_ApplyRule_NodeLabelPrefix.icc"

#endif /*GGL_MR_APPLYRULE_NODELABELPREFIX_HH_*/
