#ifndef GGL_CHEM_SMILESWRITER_HH_
#define GGL_CHEM_SMILESWRITER_HH_

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeUtil.hh"

#include <map>

namespace ggl {
 namespace chem {



 	  /*! @brief Molecule to SMILES writer
 	   *
 	   * Utility class to generate a canonical SMILES string from a molecule
 	   * graph representation. It expects atom (node) and edge (bond) label
 	   * following the Daylight's SMILES description. See 
 	   * ggl::chem::SMILES_grammar for further details.
 	   * 
 	   * It implements the algorithm suggested by Weininger (1989)
 	   *
 	   * @article{ Weininger_1989,
 	   *  title={SMILES. 2. Algorithm for generation of unique SMILES notation},
 	   *  author={Weininger, D and Weininger, A and Weininger, J L},
 	   *  journal={Journal of Chemical Information and Modeling},
 	   *  volume={29},
 	   *  number={2},
 	   *  pages={97--101},
 	   *  year={1989},
 	   *  publisher={American Chemical Society},
 	   *  url={http://dx.doi.org/10.1021/ci00062a008}
 	   * }
 	   *
 	   * NOTE : THIS WRITER IS INCOMPLETE, I.E. NOT ALL TYPES OF MOLECULE ATOMS
 	   * AND BONDS ARE HANDLED BY THIS WRITER !!!
 	   * 
 	   * Supported atom labels are defined by MoleculeUtil::getAtomData().
 	   * 
 	   * Supported bond labels are defined by MoleculeUtil::getBondData().
 	   * 
 	   * @author Alexander Ullrich (c) 2008
 	   * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
 	   * 
 	   */
 	class SMILESwriter
 	{
 	protected:
 		
 		// vertex and edge descriptor types
 		typedef boost::graph_traits<Molecule>::vertex_descriptor MVertex_t;
 		typedef boost::graph_traits<Molecule>::edge_descriptor   MEdge_t;
 		// vertex and edge iterator types
 		typedef boost::graph_traits<Molecule>::vertex_iterator     MV_Iterator_t;
 		typedef boost::graph_traits<Molecule>::adjacency_iterator  MA_Iterator_t;

 		// property map types
 		typedef boost::property_map<Molecule, PropNodeIndex>::const_type MV_Index_Map_t;
 		typedef boost::property_map<Molecule, PropNodeLabel>::const_type MV_Property_Map_t;
 		typedef boost::property_map<Molecule, PropEdgeLabel>::const_type ME_Property_Map_t;

 		//! utility map type for DFS traversal of the Molecule graph
 		typedef std::map<MVertex_t, bool> MV_Visited_Map_t;

		 /*! Container that holds all data necessary to calculate the invariant
		  * of a molecule atom node.
		  */
		class NodeData {
		public:
			MoleculeUtil::OneByte valence;
			MoleculeUtil::OneByte atomicNumber;
			MoleculeUtil::OneByte connect;
			MoleculeUtil::OneByte nonhydro;
			MoleculeUtil::OneByte sign;
			MoleculeUtil::OneByte charge;
			MoleculeUtil::OneByte protons;
			MoleculeUtil::OneByte isAromatic;
			int classLabel;
			  //! default construction
			NodeData()
			 : valence(0), atomicNumber(0), connect(0), nonhydro(0), sign(0), charge(0), protons(0), isAromatic(0), classLabel(-1)
			{}
			  /*! construction from atom label
			   *  @param atomLabel the atom label that defines the initial data
			   *  @param allowWildcard whether or not a wildcard is an allowed
			   *         atom label
			   *  @throws std::runtime_error if no atom data can be accessed for the atom label
			   */
			NodeData( const std::string & atomLabel
					, const bool allowWildcard );

			  //! calculates the node invariant based on the current data
			  //! @return the invariant of this molecule atom
			int
			getInvariant() const {
				return	  10000000*connect
						+ 100000*nonhydro
						+ 1000*atomicNumber
						+ 100*sign
						+ 10*charge
					    + protons;
			}

		};


		//! mapping of nodes to according invariant information for canonization
		typedef std::map<MVertex_t, NodeData> MV_NodeData_Map;

		static std::set<std::string> organic_subset;
 		
 	public:
 		
 		SMILESwriter();
 		
 		  /*! Generates a canonical SMILES string of the given graph 
 		   * representation of a molecule.
 		   * 
 		   * NOTE: THE FUNCTIONALITY IS INCOMPLETE, i.e. not all atom and
 		   * bond types are possible! See class description!
 		   * 
 		   * @param m the molecule graph to parse
 		   * @param ignoreProtons if true and @a m is no HH molecule,
 		   *        all protons that can be inferred from atom and bond valence
 		   *        are ignored when producing the SMILES string;
 		   *        otherwise protons are compressed
 		   *        into the adjacent non-proton atom label
 		   * @param allowWildcard whether or not the wildcard is a valid atom
 		   *        label within SMILES. NOTE: this will result in non-standard
 		   *        SMILES since the wildcard label is a non-standard extension
 		   *        within the GGL!
 		   * @return a canonical SMILES string representing the given molecule
		   * @throws std::runtime_error if unsupported atom or bond labels are
		   *        encountered
 		   */
 		static
 		std::string
 		getSMILES( const Molecule& m
 					, const bool ignoreProtons = true
 					, const bool allowWildcard = false );
 		
 		  /*! Generates a canonical SMILES string of the given graph
 		   * representation of a molecule.
 		   *
 		   * Before this is done, all known groups are compressed into their
 		   * according group labels using MoleculeUtil::compressGroups(..).
 		   *
 		   * NOTE: THE FUNCTIONALITY IS INCOMPLETE, i.e. not all atom and
 		   * bond types are possible! See class description!
 		   *
 		   * @param m the molecule graph to parse
	       * @param groups a container that holds group IDs where
	       *        each matching node represents the according subgraph
 		   * @param ignoreProtons if true and @a m is no HH molecule,
 		   *        all protons that can be inferred from atom and bond valence
 		   *        are ignored when producing the SMILES string;
 		   *        otherwise protons are compressed
 		   *        into the adjacent non-proton atom label
 		   * @param allowWildcard whether or not the wildcard is a valid atom
 		   *        label within SMILES. NOTE: this will result in non-standard
 		   *        SMILES since the wildcard label is a non-standard extension
 		   *        within the GGL!
 		   * @return a canonical SMILES string representing the given molecule
 		   */
 		static
 		std::string
 		getSMILES( const Molecule& m
 				, const GroupMap & groups
 				, const bool ignoreProtons
				, const bool allowWildcard = false );

 	protected:
 		
 		  /*! Generates a canonical SMILES string of the given graph
 		   * representation of a molecule.
 		   *
 		   * NOTE: THE FUNCTIONALITY IS INCOMPLETE, i.e. not all atom and
 		   * bond types are possible! See class description!
 		   *
 		   * @param m the molecule graph to parse
	       * @param groups a container that holds group IDs where
	       *        each matching node represents the according subgraph; can be
	       *        NULL if no groups are to be considered
 		   * @param ignoreProtons if true and @a m is no HH molecule,
 		   *        all protons that can be inferred from atom and bond valence
 		   *        are ignored when producing the SMILES string;
 		   *        otherwise protons are compressed
 		   *        into the adjacent non-proton atom label
 		   * @param allowWildcard whether or not the wildcard is a valid atom
 		   *        label within SMILES. NOTE: this will result in non-standard
 		   *        SMILES since the wildcard label is a non-standard extension
 		   *        within the GGL!
 		   * @return a canonical SMILES string representing the given molecule
 		   */
 		static
 		std::string
 		getSMILES( const Molecule& m
 				, const GroupMap * const groups
 				, const bool ignoreProtons
				, const bool allowWildcard );

 		static
 		const int primes[];
 		
 		static
 		const int primesLength;
 		
 		static 
 		int 
 		prime(int number);
 		
 		static
 		std::set<std::string>
 		initOrganicSubset();

 		 /*
	      * @param groups a container that holds group IDs where
	      *        each matching node represents the according subgraph; can be
	      *        NULL if no groups are to be considered
	      */
 		static 
 		std::vector<int>
 		canonize(const Molecule *graph
 				, const MV_Visited_Map_t &visit
				, const MV_Index_Map_t & idx
				, const MV_Property_Map_t & vname
				, const ME_Property_Map_t & ename
				, MV_NodeData_Map & nodeData
 				, const GroupMap * const groups
 				, const bool allowWildcard
 				);
 		
 		static
 		std::string
 		build_dfs(	MVertex_t v, MVertex_t p
 					, const Molecule *graph
 					, std::vector<int> *ranks
 					, MV_Visited_Map_t & visit
					, const MV_NodeData_Map & nodeData
					, const MV_Index_Map_t & idx
					, const MV_Property_Map_t & vname
					, const ME_Property_Map_t & ename
					, const bool ignoreProtons
 					);

 		static
 		std::string 
 		second_pass(std::string smiles);
 		
		  /*!
		   * Checks whether or not an atom identifier has to be enclosed in
		   * brackets within the SMILES notation or not.
		   *
		   * @param atom the atom label without enclosing brackets
		   * @return true, if the atom label has to be enclosed with brackets;
		   *         false otherwise.
		   */
		static
		bool
		isWithBracketsInSMILES( const std::string& atom );

		  /*!
		   * Produces a SMILES conform label of the given atom label, ie.
		   * enclosing brackets are added if needed.
		   * @param label the atom label to check
		   * @return the SMILES conform atom label
		   */
		static
		std::string
		getLabel( const std::string& label );

		  /*!
		   * Produces a SMILES conform label of the given atom label, ie.
		   * enclosing brackets are added if needed. Furthermore, protons can
		   * be removed from the SMILES if their number can be deduced from
		   * atom and bond valence.
		   * @param label the atom label to check
		   * @param nodeData the atom information for this atom
		   * @param ignoreProtons whether or not protons should be removed from
		   *        the atom label if possible
		   * @return the SMILES conform atom label
		   */
		static
		std::string
		getLabel( const std::string& label
				, const NodeData& nodeData
				, const bool ignoreProtons );

 	};
 	
 
 } // namespace chem
} // namespace ggl

  // include method implementation
#include "ggl/chem/SMILESwriter.icc"

#endif /*SMILESWRITER_HH_*/
