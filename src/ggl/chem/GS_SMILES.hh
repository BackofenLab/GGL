#ifndef GGL_CHEM_GS_SMILES_HH_
#define GGL_CHEM_GS_SMILES_HH_

#include "ggl/chem/GS_chem.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESwriter.hh"

#include <iterator>

namespace ggl {
 namespace chem {

 
#define GGL_CHEM_GS_SMILES_TEMPLATE \
	template <	class STL_INSERTER >

#define GGL_CHEM_GS_SMILES_TYPE \
 	GS_SMILES <	STL_INSERTER >

 
	  /*! @brief SMILES graph storage
	   *
	   * A Graph_Storage implementation that converts each added Molecule graph
	   * into a SMILES string representation and adds it to the specified 
	   * STL string container.
	   * 
	   * Usage Examples :
	   * \verbatim
	   * //////////////////////////////////////////////////////////////////////
	   * 
	   * std::set< std::string > SMILES_set;
	   * typedef std::insert_iterator< SMILES_set > SMILES_set_inserter;
	   * 
	   * SMILES_set smilesSet;
	   * SMILES_set_inserter insertSet(smilesSet, smilesSet.end());
	   * 
	   * ggl::chem::GS_SMILES< SMILES_set_inserter > gs_set(insertSet);
	   * 
	   * //////////////////////////////////////////////////////////////////////
	   * 
	   * std::vector< std::string > SMILES_vector;
	   * typedef std::back_insert_iterator< SMILES_vector > SMILES_vector_inserter;
	   * 
	   * SMILES_vector smilesVec;
	   * SMILES_vector_inserter insertVec(smilesVec);
	   * 
	   * ggl::chem::GS_SMILES< SMILES_vector_inserter > gs_vec(insertVec);
	   * 
	   * //////////////////////////////////////////////////////////////////////
	   * \endverbatim
	   * 
	   * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   * 
	   * @tparam STL_INSERTER an STL insert iterator (e.g. std::insert_iterator)
	   *            to add all SMILES to store to
	   * @tparam Molecule the boost graph representing the molecule to convert
	   * @tparam NODE_LABEL_PROPERTY property type of the node label
	   * @tparam EDGE_LABEL_PROPERTY property type of the edge label
	   * @tparam NODE_INDEX_PROPERTY property type of the node indices
	   */
	template <	class STL_INSERTER >
	class GS_SMILES : public GS_chem {
		
	protected:
		
		  //! the std::inserter where each SMILES string is reported to
		STL_INSERTER insert;
		

	public:
		
		  //! Construction
		  //! @param insert the STL inserter to which each generated SMILES
		  //!               is assigned to
		GS_SMILES( STL_INSERTER insert );
		
		
		virtual
		~GS_SMILES();
		
		
		  //! Writes the SMILES string of a given molecule graph to a string
		  //! container.
		  //! @param graph the molecule object to add.
		virtual
		void
		addMolecule( const Molecule & graph );
	};

 
 }  // namespace chem
}  // namespace ggl


namespace ggl {
 namespace chem {
 
#define GGL_CHEM_GS_SMILES_MOL_TEMPLATE \
 	template <	class SMILES_MOL_MAP >

#define GGL_CHEM_GS_SMILES_MOL_TYPE \
 	GS_SMILES_MOL <	SMILES_MOL_MAP >


	 /*! @brief SMILES graph storage using STL SMILES-Graph map
	  *
	  * A Graph_Storage implementation that converts each added Molecule graph
	  * into a SMILES string representation and adds it, if not already
	  * existing, to the specified STL map container using the SMILES as key
	  * and the Molecule object as value.
	  *
	  * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  * @tparam STL_INSERTER an STL insert iterator (e.g. std::insert_iterator)
	  *            to add all SMILES to store to
	  */
	template <	class SMILES_MOL_MAP >
	class GS_SMILES_MOL : public GS_chem {
		
	protected:
		
		  //! the map where each SMILES string is mapped to the represented
		  //! molecule
		SMILES_MOL_MAP* smiles2mol;

	public:
		
		  //! Construction
		  //! @param smiles2mol the STL inserter to which each generated SMILES
		  //!               and its molecule is assigned to
		GS_SMILES_MOL( SMILES_MOL_MAP& smiles2mol );
		
		
		virtual
		~GS_SMILES_MOL();
		
		
		  //! Converts a given molecule graph to SMILES and adds it to the
		  //! storage container.
		  //! @param graph the Graph object to add.
		virtual
		void
		addMolecule( const Molecule & graph );
		
	protected:
		
		virtual
		bool
		insert2map(const std::string & SMILES, const Molecule& graph );
		
	};

 }  // namespace chem
}  // namespace ggl


namespace ggl {
 namespace chem {

#define GGL_CHEM_GS_SMILES_MOLp_TEMPLATE \
 	template <	class SMILES_MOL_MAP >

#define GGL_CHEM_GS_SMILES_MOLp_TYPE \
 	GS_SMILES_MOLp <	SMILES_MOL_MAP >

	 /*! @brief SMILES graph storage using STL SMILES-GraphPointer map
	  *
	  * A Graph_Storage implementation that converts each added Molecule graph
	  * into a SMILES string representation and adds it, if not already
	  * existing, to the specified STL map container using the SMILES as key
	  * and the pointer to the newly allocated Molecule object as value.
	  *
	  * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  * @tparam STL_INSERTER an STL insert iterator (e.g. std::insert_iterator)
	  *            to add all SMILES to store to
	  */
	template <	class SMILES_MOL_MAP >
	class GS_SMILES_MOLp : public GS_chem {
		
	protected:
		
		  //! the map where each SMILES string is mapped to the represented
		  //! molecule
		SMILES_MOL_MAP* smiles2mol;
		const SMILES_MOL_MAP* smiles2molCheckOnly1;
		const SMILES_MOL_MAP* smiles2molCheckOnly2;


	public:
		
		  //! Construction
		  //! @param smiles2mol the map (SMILES->molecule) to that each 
		  //!               generated SMILES and its molecule is added to
		  //!               if not already present
		GS_SMILES_MOLp( SMILES_MOL_MAP& smiles2mol );
		
		  /*! Construction
		   * @param smiles2mol the map (SMILES->molecule) to that each
		   *                generated SMILES and its molecule is added to
		   *                if not already present
		   * @param smiles2molCheckOnly an additional map to check for existence
		   *                of a new molecule to store
		   */ 
		GS_SMILES_MOLp( SMILES_MOL_MAP& smiles2mol 
						, const SMILES_MOL_MAP& smiles2molCheckOnly );
		
		  /*! Construction
		   * @param smiles2mol the map (SMILES->molecule) to that each
		   *                generated SMILES and its molecule is added to
		   *                if not already present
		   * @param smiles2molCheckOnly1 an additional map to check for existence
		   *                of a new molecule to store
		   * @param smiles2molCheckOnly2 an additional map to check for existence
		   *                of a new molecule to store
		   */ 
		GS_SMILES_MOLp( SMILES_MOL_MAP& smiles2mol 
						, const SMILES_MOL_MAP& smiles2molCheckOnly1
						, const SMILES_MOL_MAP& smiles2molCheckOnly2 );
		
		
		virtual
		~GS_SMILES_MOLp();
		
		
		  //! Converts a given molecule graph to SMILES and adds it to the
		  //! storage container.
		  //! @param graph the Graph object to add.
		virtual
		void
		addMolecule( const Molecule & graph );
		
	protected:
		
		virtual
		bool
		insert2map(const std::string & SMILES, const Molecule& graph );
		
	};

 }  // namespace chem
}  // namespace ggl

#include "ggl/chem/GS_SMILES.icc"


#endif /*GGL_CHEM_GS_SMILES_HH_*/
