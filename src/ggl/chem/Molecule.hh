#ifndef GGL_CHEM_MOLECULE_HH_
#define GGL_CHEM_MOLECULE_HH_

#include "ggl/Graph.hh"
#include "sgm/Graph_boost.hh"
#include "sgm/Graph_boostV_p.hh"
	
namespace ggl {
  namespace chem {


	  //! This boost graph property is used to determine the index of a given
	  //! atom node along the iterator order.
	typedef ggl::PropNodeIndex	PropNodeIndex;

	  //! Vector of node indices
	typedef ggl::NodeIndexVec NodeIndexVec;

	  //! Set of node indices
	typedef ggl::NodeIndexSet NodeIndexSet;

	  //! This boost graph property is used to determine the label of a given
	  //! atom node.
	typedef ggl::PropNodeLabel	PropNodeLabel;

	  //! This boost graph property is used to determine the label of a given
	  //! bond edge.
	typedef ggl::PropEdgeLabel	PropEdgeLabel;


	  //! The boost properties stored for an atom in a molecule
	  //!
	  //! @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	typedef ggl::Graph_NodeProperties Molecule_AtomProperties;
	
	  //! The boost properties stored for a bond in a molecule
	  //!
	  //! @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	typedef ggl::Graph_EdgeProperties Molecule_BondProperties;
	
	  //! The boost graph based molecule representation of atoms and bonds
	  //!
	  //! @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	typedef ggl::Graph Molecule;


	  //! Wrapper typedef to use a Molecule within the SGM Graph_Interface
	typedef sgm::Graph_boost< Molecule
							, PropNodeLabel
							, PropEdgeLabel
							, PropNodeIndex >	Molecule_Graph;


	  //! Wrapper typedef to use a vector of Molecule objects within the
	  //! SGM Graph_Interface
	typedef sgm::Graph_boostV_p< Molecule
							, PropNodeLabel
							, PropEdgeLabel
							, PropNodeIndex >	Molecule_Graph_V;
  
  }
}

#endif /*MOLECULE_HH_*/
