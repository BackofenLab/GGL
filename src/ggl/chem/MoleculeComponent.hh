
#ifndef GGL_CHEM_MOLECULECOMPONENT_HH_
#define GGL_CHEM_MOLECULECOMPONENT_HH_

#include <string>
#include <vector>
#include <vector>

#include "sgm/Pattern.hh"
#include "ggl/Graph.hh"

namespace ggl {
namespace chem {

	/*! @brief Molecule component definition
	 *
	 * Describes a molecule component that e.g. is used by
	 * ggl::chem::MoleculeDecomposition
	 * to assess thermodynamic properties of a molecule. It describes the
	 * molecule component pattern to be matched and the according
	 * properties.
	 *
	 * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	class MoleculeComponent {

	public:

		  //! type of the pattern graph
		typedef ggl::Graph PatternGraph;

		  //! type of a node ID set
		typedef std::set<size_t> NodeSet;

		  //! type of a Match_Constraint vector
		typedef sgm::Pattern_Interface::ConstraintVec ConstraintVec;

		  //! type of a ring fragment definition, i.e. a list of nodes
		typedef std::vector< size_t > RingFragmentList;

		  //! different types of ring fragments to consider
		enum RingFragmentType {
			RF_aromaticHydrocarbon = 0
			, RF_heteroaromatic = 1
			, RF_nonaromatic = 2
			, RF_undefined = 3
		};

		  //! a ring fragment description
		class RingFragment {
		public:
			  //! the type of the ring fragment
			RingFragmentType type;
			  //! the definition of the ring fragment node list
			RingFragmentList fragment;

			  //! construction
			  //! @param type the type to set
			  //! @param fragment the fragment definition to set
			RingFragment( const RingFragmentType& type
						, const RingFragmentList& fragment);
		};

		  //! container for ring fragments
		typedef std::vector< RingFragment > RingFragVec;

	public:
		  //! description of the component to be matched
		std::string description;
		  //! the priority of this component
		  //! (used in ggl::chem::MoleculeDecomposition for an ordering)
		size_t priority;
		  //! the free energy contribution of the component
		  //! (used in ggl::chem::MoleculeDecomposition for energy calculations)
		double freeEnergy;
		  //! the pattern to be matched for this component
		PatternGraph pattern;
		  //! the set of node IDs from the pattern that are described by this
		  //! component
		NodeSet compIDs;
		  //! additional matching constraints
		ConstraintVec constraints;
		  //! additional ring fragments information that has to be matched
		RingFragVec ringFragments;

	public:
		  //! default construction
		MoleculeComponent();

		  //! copy construction
		  //! @param toCopy the object to make this a copy of
		MoleculeComponent( const MoleculeComponent& toCopy );

		  //! destruction
		~MoleculeComponent();

		  //! assignment operator
		  //! @param toCopy the object to make this a copy of
		  //! @return access to the changed *this object
		MoleculeComponent&
		operator = (const MoleculeComponent& toCopy );

	};

}  // namespace chem
} // namespace ggl

#endif /* GGL_CHEM_MOLECULECOMPONENT_HH_ */
