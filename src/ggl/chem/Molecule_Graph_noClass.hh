#ifndef GGL_CHEM_MOLECULE_GRAPH_NOCLASS_HH_
#define GGL_CHEM_MOLECULE_GRAPH_NOCLASS_HH_

#include "sgm/Graph_NodeLabelPrefix.hh"

namespace ggl {
  namespace chem {

	  /*! @brief Wrapper for a molecule representing graph without atom class label
	   *
	   *  Represents a molecule class where class information is removed from
	   *  all atom labels
	   *
	   *  @author Martin Mann (c) 2017 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Molecule_Graph_noClass : public sgm::Graph_NodeLabelPrefix {
		
	public:

		//! Construction
		//! @param molecules the graph representing the molecules to wrap;
		//!    atom labels are reduced to their substrings up to the first
		//!    occurrence of the class label separator ':'
		Molecule_Graph_noClass( const sgm::Graph_Interface & molecules )
		 :	Graph_NodeLabelPrefix(molecules, ":")
		{}

		  //! Destruction
		virtual 
		~Molecule_Graph_noClass() {}

	};

  } // namespace chem
} // namespace ggl



#endif /*GGL_CHEM_MOLECULE_GRAPH_NOCLASS_HH_*/
