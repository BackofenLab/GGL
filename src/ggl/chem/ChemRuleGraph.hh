#ifndef GGL_CHEM_CHEMRULEGRAPH_HH_
#define GGL_CHEM_CHEMRULEGRAPH_HH_

#include "ggl/chem/ChemRule.hh"
#include "ggl/RuleGraph.hh"

namespace ggl {
 namespace chem {
 
 
 	
 	  //! RightSidePattern graph of a ggl::chem::ChemRule
 	typedef ggl::RightSidePattern RightSidePattern;
 	
 	
 	  //! LeftSidePattern graph of a ggl::chem::ChemRule
 	class LeftSidePattern : public ggl::LeftSidePattern
 	{

 	public:

 		  //! type of the super class
 		typedef ggl::LeftSidePattern SuperClass;

 		  //! index type used
		typedef SuperClass::IndexSet IndexSet;

		  //! Constructs a LeftSidePattern for the given rule
		  //! @param rule the ChemRule this object will be a left side pattern of
		LeftSidePattern( const ChemRule& rule )
		 : SuperClass(rule)
		{}

		  //! destruction
		virtual
		~LeftSidePattern( )
		{}

		  //! Constant access to the internal ChemRule object.
		  //! @return the internal Rule object
		virtual
		const ChemRule &
		getRule(void) const {
			return (const ChemRule&)(SuperClass::getRule());
		}


 	};

 } // chem
} // ggl

#endif /*GGL_CHEM_CHEMRULEGRAPH_HH_*/
