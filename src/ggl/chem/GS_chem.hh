#ifndef GGL_CHEM_GS_CHEM_HH_
#define GGL_CHEM_GS_CHEM_HH_

#include "ggl/Graph_Storage.hh"

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>

#include "ggl/chem/Molecule.hh"

namespace ggl {
 namespace chem {




	/*! @brief Interface molecule graph storage
	 *
	 * An abstract interface for the storage of molecule objects.
	 * For instance, it is used by ggl::MR_ApplyRule to store the
	 * molecule graphs resulting from graph grammar ggl::chem::ChemRule
	 * applications.
	 *
	 * @author Martin Mann (c) 2011 http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	class GS_chem : public ggl::Graph_Storage
	{
	protected:

		  //! the container type used to represent the connected component id
		  //! for each node
		typedef std::vector<int> ComponentIdVec;

		  //! the type of a constant index map of a Molecule boost graph
		typedef boost::property_map<Molecule, PropNodeIndex>::const_type IndexMap;

		  /*!
		   * A boost graph node filter used as VertexPredicate of a
		   * boost::filtered_graph to represent a specific connected component
		   * of a Molecule. This is needed if the added Molecule represents
		   * a set of molecules to generate a SMILES string for each of them.
		   *
		   */
	 	class node_is_in_component {
	 	protected:
	 		  //! componentID the index of the connected component to represent
	 		int componentID;
	 		  //! the connected component ID for each vertex
	 		const ComponentIdVec* componentVec;
	 		  //! the index map of the underlying graph to get the index of a
	 		  //! requested vertex
	 		const IndexMap* idxMap;
	 	public:
	 		  //! empty construction
	 		node_is_in_component();

	 		  /*! construction
	 		   * @param componentID the index of the connected component to
	 		   *        represent
	 		   * @param componentVec the connected component ID for each vertex
	 		   * @param idxMap the index map of the underlying graph to get the
	 		   *        index of a requested vertex
	 		   */
	 		node_is_in_component(	const int componentID
	 								, const ComponentIdVec& componentVec
	 								, const IndexMap* idxMap );

	 		  //! Check if a vertex should be part of the filtered graph that
	 		  //! represents a specific connected component of a boost graph.
	 		  //! @return true if the node should be represented, false otherwise.
	 		template <typename VERTEX>
	 		bool operator()(const VERTEX& node) const;
	 	};

		  /*!
		   * A boost graph edge filter used as EdgePredicate of a
		   * boost::filtered_graph to represent a specific connected component
		   * of a Molecule. This is needed if the added Molecule represents
		   * a set of molecules to generate a SMILES string for each of them.
		   *
		   */
	 	class edge_is_in_component {
	 	protected:
	 		  //! componentID the index of the connected component to represent
	 		int componentID;
	 		  //! the connected component ID for each vertex
	 		const ComponentIdVec* componentVec;
	 		  //! the index map of the underlying graph to get the index of a
	 		  //! requested vertex
	 		const IndexMap* idxMap;
	 		  //! the molecule graph where the edges of interest are from
	 		const Molecule* m;
	 	public:

	 		 //! empty construction
	 		edge_is_in_component();

	 		  /*! construction
	 		   * @param componentID the index of the connected component to
	 		   *        represent
	 		   * @param componentVec the connected component ID for each vertex
	 		   * @param idxMap the index map of the underlying graph to get the
	 		   *        index of a requested vertex
	 		   * @param m the molecule graph where the edges of interest are
	 		   *        from
	 		   */
	 		edge_is_in_component(	const int componentID
	 								, const ComponentIdVec& componentVec
	 								, const IndexMap* idxMap
	 								, const Molecule& m);


	 		  //! Check if an edge should be part of the filtered graph that
	 		  //! represents a specific connected component of a boost graph.
	 		  //! @return true if source and target graph are within the
	 		  //!        component and should be represented, false otherwise.
	 		template <typename EDGE>
	 		bool operator()(const EDGE& e) const;
	 	};


		  //! a component graph definition if more than one connected
	 	  //! component present is present and it is needed to split the graph
	 	  //! to report into its components
		typedef boost::filtered_graph<	Molecule
												, edge_is_in_component
												, node_is_in_component
											> ComponentGraph;

	public:
		
		virtual ~GS_chem()
		{}
		
		  //! The reported graphs is split into its individual independent
		  //! components. Each component is forwarded to addMolecule to be
		  //! implemented by any sub class.
		  //! @param graph the graph object to add that encodes one or more
		  //!        molecules
		virtual
		void
		add( const Molecule & graph );

		  //! Adds a single molecule graph to the storage.
		  //! @param mol the molecule graph to add
		virtual
		void
		addMolecule( const Molecule & mol ) = 0;

	};

 } // namespace chem
} // namespace ggl

#include "ggl/chem/GS_chem.icc"

#endif /*GGL_CHEM_GS_CHEM_HH_*/

