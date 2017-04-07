#ifndef GGL_CHEM_MOLECULEUTIL_HH_
#define GGL_CHEM_MOLECULEUTIL_HH_

#include <map>
#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <sgm/HashMap.hh>

#if HAVE_UNORDERED_MAP > 0
	#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP > 0
	#include <tr1/unordered_map>
#elif HAVE_GNU_HASH_MAP > 0
	#include <ext/hash_map>
#else
	#include <map>
#endif

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeComponent.hh"

namespace ggl {
	namespace chem {

	
	 /*! Container for molecule group ID -> subgraph mappings
	  *
	  * Container to store the IDs of molecular groups that can be represented
	  * within molecules by a single node and that are replaced by according
	  * subgraphs (stored within this container).
	  */
	typedef
	#if HAVE_UNORDERED_MAP > 0
		std::unordered_map<std::string, MoleculeComponent >
	#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map<std::string, MoleculeComponent >
	#elif HAVE_GNU_HASH_MAP > 0
		__gnu_cxx::hash_map< std::string, MoleculeComponent, sgm::hash_string >
	#else
		std::map< std::string, MoleculeComponent >
	#endif
		GroupMap;

	  /*! @brief Molecule utility factory
	   *
	   * Utility class that contains certain data and members needed for
	   * Molecule evaluation and handling. 
	   * 
 	   * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class MoleculeUtil
	{
	public:
		
		  //! typedef to represent one byte of information
 		typedef unsigned char OneByte;

		  //! the wildcard character for atom labels valid in molecule
		  //! descriptions
		static const std::string AtomLabelWildcard;
		
		 /*! Data type that stores information connected to a certain atom label.
		  */
		class AtomLabelData {

		protected:

			static std::vector<OneByte> shellSizeSums;

		public:

			typedef std::vector<OneByte> OneByteVec;

			  //! atomic number
			OneByte atomicNumber;
			  //! available standard number of valence electrons
			OneByte valence;
			  //! is set to 1 if the label denotes a part of an aromatic ring, 
			  //! 0 otherwise
			OneByte isAromatic;
			  //! mean atomic weight
			double atomicWeight;
			  //! maximal number of protons attached
			OneByte maxProtons;
			  //! list of all possible numbers of "valences" to be considered for valence checks
			OneByteVec valenceAlternatives;
			  //! if set to 1 (default) a valence check is done for this atom,
			  //! otherwise (0) no valence check is performed
			OneByte isToBeChecked;
		  
			  /*! Construction
			   * initializes the list of all possible valences with "valence"
			   *
			   * @param atomicNumber_ the atomic number of the atom
			   * @param valence_ the available number of valence electrons
			   * @param isAromatic_ set 1 if the atom is part of an aromatic
			   *         ring; 0 otherwise
			   * @param atomicWeight_ the atomic weight of the atom
			   * @param isToBeChecked set 1 if a valence check is to be done
			   *         for this type of atom, 0 otherwise.
			   */
			AtomLabelData(	const OneByte atomicNumber_
							, const OneByte valence_
							, const OneByte isAromatic_
							, const double atomicWeight_
							, const OneByte isToBeChecked_ = 1)
			 :	  atomicNumber(atomicNumber_)
				, valence(valence_)
				, isAromatic(isAromatic_)
			 	, atomicWeight(atomicWeight_)
				, maxProtons(valence)
				, valenceAlternatives(1)
				, isToBeChecked(isToBeChecked_)
		 	{
				  // initialize list of all valences to be considered for checks
				*(valenceAlternatives.begin()) = valence;
		 	}

			  /*! Construction
			   * @param atomicNumber_ the atomic number of the atom
			   * @param valence_ the number of valence electrons to set
			   * @param isAromatic_ set 1 if the atom is part of an aromatic
			   *         ring; 0 otherwise
			   * @param atomicWeight_ the atomic weight of the atom
			   * @param valenceList list of all valences to be considered for checks.
			   * @param isToBeChecked set 1 if a valence check is to be done
			   *         for this type of atom, 0 otherwise.
			   */
			AtomLabelData(
							const OneByte atomicNumber_
							, const OneByte valence_
							, const OneByte isAromatic_
							, const double atomicWeight_
							, const OneByte isToBeChecked_
							, const OneByteVec valenceList )
			 :	  atomicNumber(atomicNumber_)
				, valence(valence_)
				, isAromatic(isAromatic_)
			 	, atomicWeight(atomicWeight_)
				, valenceAlternatives(valenceList)
				, isToBeChecked(isToBeChecked_)
		 	{
				  // ensure that the standard valence is among the possible valence values
				if( std::find(valenceAlternatives.begin(), valenceAlternatives.end(), valence) == valenceAlternatives.end() ) {
					valenceAlternatives.push_back(valence);
				}
				  // does not make sense to give alternatives and don't enable checking ...
				assert( valenceAlternatives.size() == 1 || (int)isToBeChecked == 1 );
		 	}

		};
		
		  /*! Mapping of atom labels (in SMILES notation but without brackets)
		   * to the corresponding atom information like valence etc.
		   */
		typedef std::map< std::string, AtomLabelData > AtomDataMap;

		  /*! Mapping of aromatic labels to their non-aromatic form and vice
		   * versa, to enable a relabeling of aromatic rings.
		   */
		typedef std::map< std::string, std::string > AromaticSwapMap;
		
		 /*! Data type that stores information connected to a certain bond label.
		  */
		class BondLabelData {
		public:
			  //! number of valence electrons of the bond
			OneByte valence;
			  //! is set to 1 if the label denotes a part of an aromatic ring, 
			  //! 0 otherwise
			OneByte isAromatic;
			
			  /*! Construction
			   * @param valence_ the number of valence electrons to set
			   * @param isAromatic_ set 1 if the bond is part of an aromatic ring;
			   *         0 otherwise
			   */
			BondLabelData(	const OneByte valence_
							, const OneByte isAromatic_ )
			 : valence(valence_)
			 	, isAromatic(isAromatic_)
		 	{}
		};

		  /*! Mapping of bond labels (in SMILES notation) to the corresponding
		   * bond information like valence etc.
		   */
		typedef std::map< std::string, BondLabelData > BondDataMap;


	public:

		  //! consistency code : everything fine 
		static const size_t C_Consistent;
		  //! consistency code : at least one atom label is not SMILES conform
		  //! or currently not supported within the library
		static const size_t C_AtomLabelInvalid;
		  //! consistency code : at least one bond label is not SMILES conform
		  //! or currently not supported within the library
		static const size_t C_BondLabelInvalid;
		  //! consistency code : at least one atom label is a wildcard which is
		  //! currently not supported within the library
		static const size_t C_AtomLabelWildcard;
		  //! consistency code : at least one atom label is complex and contains
		  //! implicit H atoms, this is currently not supported
		static const size_t C_AtomComplexWithH;
		  //! consistency code : at least one atom shows an inconsistent
		  //! electron distribution, ie. valence+charge != protonCount
		static const size_t C_AtomValence;
		  //! consistency code : at least one bond connects the same atom, i.e.
		  //! forms a loop
		static const size_t C_BondLoop;
		  //! consistency code : the molecule graph is not connected, i.e.
		  //! contains at least two connected components
		static const size_t C_NonConnected;

	public:
		
		  //! Default construction
		MoleculeUtil();
		  //! Default destruction
		virtual ~MoleculeUtil();
		
 		///////////////////  ATOM DATA ACCESS ETC.  ////////////////////////
		
		  /*! Access to the currently supported atom labels and the 
		   * corresponding atom information, e.g. valence etc.
		   * 
		   * @return the atom2data mapping
		   */
 		static
 		const AtomDataMap&
 		getAtomData( void );
		
		  /*! Access to atom information for the given atom label, if existing.
		   * 
		   * @param label the atom label to derive the information from
		   * @return the atom2data entry or NULL if none available
		   */
 		static
 		const AtomLabelData * const
 		getAtomData( const std::string& label );
 		

 		  /*!
 		   * Access to the aromatic/non-aromatic pendant of an atom or edge
 		   * label if it exists.
 		   *
 		   * @param label the label of interest
 		   * @return a pointer to the pendant label or NULL if none exists
 		   */
 		static
 		const std::string * const
 		getAromaticPendant( const std::string& label );

		  /*! Access to atom name within a given atom label.
		   *
		   * @param label the atom label to derive the information from
		   * @return the name of the atom
		   */
 		static
 		std::string
 		getAtom( const std::string& label );

		  /*! Access to number of protons within a given atom label,
		   * if existing.
		   *
		   * @param label the atom label to derive the information from
		   * @return the number of additional protons in the label
		   */
 		static
 		size_t
 		getProtons( const std::string& label );
 		
		  /*! Access to the charge within a given atom label,
		   * if existing.
		   *
		   * @param label the atom label to derive the information from
		   * @return the charge information in the label
		   */
 		static
 		int
 		getCharge( const std::string& label );
 		
		  /*! Access to the class information within a given atom label,
		   * if existing.
		   *
		   * @param label the atom label to derive the information from
		   * @return the class information in the label, or 0 if not present
		   */
 		static
 		int
 		getClass( const std::string& label );
 		

 		  /*!
 		   * Produces a complex atom label with the given information.
 		   *
 		   * @param atom the atom label
 		   * @param protons the number of protons attached to the atom;
 		   *        a value of 0 is ignored
 		   * @param charge the charge of the atom; a value of 0 is ignored
 		   * @param classID the classID; a value of 0 is ignored
 		   * @param explicitChargeValue if true, a charge value of 1 is
 		   * 		represented by "+1" rather than just "+"
 		   *
 		   * @return the according complex atom label
 		   */
 		static
 		std::string
 		getComplexAtomLabel( const std::string& atom
							, const size_t protons = 0
							, const int charge = 0
							, const int classID = 0
							, const bool explicitChargeValue = false );


 		///////////////////  BOND DATA ACCESS ETC.  ////////////////////////
		
		  /*! Access to the currently supported bond labels and the 
		   * corresponding bond information, e.g. valence etc.
		   * 
		   * @return the bond2data mapping
		   */
 		static
 		const BondDataMap&
 		getBondData( void );
		
		  /*! Access to bond information for the given bond label, if existing.
		   * 
		   * @param label the bond label to derive the information from
		   * @return the bond2data entry or NULL if none available
		   */
 		static
 		const BondLabelData * const
 		getBondData( const std::string& label );
 		
 		///////////////////  CONSISTENCY CHECKS ETC.  ////////////////////////
		
		  /*! Checks if the given atom (node) label is SMILES conform and 
		   * currently supported, e.g. by the SMILESwriter class. 
		   * 
		   * @param atomLabel the label to check
		   * @return true if the atom label is ok, false otherwise
		   */
		static
		bool
		isValidAtomLabel( const std::string & atomLabel );
		
		  /*! Checks if the given bond (edge) label is SMILES conform and 
		   * currently supported, e.g. by the SMILESwriter class. 
		   * 
		   * @param bondLabel the label to check
		   * @return true if the bond label is ok, false otherwise
		   */
		static
		bool
		isValidBondLabel( const std::string & bondLabel );
 		
 		  /*! Checks if a given ggl::chem::Molecule graph is consistent. 
 		   * If not, according error codes are returned.
 		   * 
 		   * @param mol the molecule graph to check
 		   * @return if consistent C_Consistent is returned. Otherwise a product
 		   *        of the according C_* values.
 		   */
 		static
 		size_t
 		isConsistent( const Molecule & mol );
 		
 		  /*!
 		   * Writes a description of the consistency status or errors, encoded
 		   * in a consistency code produced by isConsistent*(...), to a given
 		   * outstream. The function returns whether or not an error occured.
 		   *
 		   * @param consistencyCode the error code to parse, produced by
 		   *           a call to isConsistent*(...)
 		   * @param errorStream the output stream to write the error description
 		   *           to
 		   * @return true if no error is encoded; false otherwise
 		   */
 		static
 		bool
 		decodeConsistencyStatus(	const size_t consistencyCode
									, std::ostream& errorStream ) ;


 		///////////////////  MOLECULE EDITING ETC.  ////////////////////////
 		
 		
 		  /*! Copies a given molecule into another molecule object which is 
 		   * overwritten.
 		   * 
 		   * @param mol the molecule graph to copy
 		   * @param toFill the molecule graph to make a copy of mol
 		   */
 		static
 		void
 		copy( const Molecule& mol, Molecule & toFill);


 		  /*! Copies a given molecule into another molecule object which is
 		   * overwritten.
 		   *
 		   * Note, no sanity checks of node/edge labels, degrees, etc. are done!
 		   *
 		   * @param mol the molecule graph to copy
 		   * @param toFill the molecule graph to make a copy of mol
 		   */
 		static
 		void
 		copy( const sgm::Graph_Interface& mol, Molecule & toFill);

 		  /*! Compresses all explicitly represented "H" atoms into the adjacent
 		   * atom label.
 		   * 
 		   * Note: only protons with atom node label "H" are compressed. Nodes
 		   * e.g. including class information like "H:1" are maintained and not
 		   * removed.
 		   *
 		   * @param mol the molecule graph to compress
 		   */
 		static
 		void
 		compressHnodes( Molecule& mol );


 		  /*! Removes all represented "H" atoms that can be inferred from
 		   * atom valence, charge and adjacent bond information. All other
 		   * proton information that cannot be inferred is collapsed into
 		   * accoring complex atom labels.
 		   *
 		   * Thus, in the end, the molecule won't show adjacent protons if any.
 		   *
 		   * Single protons or HH molecules are preserved as that.
 		   *
 		   * @param mol the molecule graph to be stripped
 		   */
 		static
 		void
 		removeProtons( Molecule& mol );


 		 /*!
 		  * Computes the number of missing hydrogens to be added to this atom.
 		  */
 		static
 		size_t
 		getProtonsToAdd( const std::string & atomLabel
 						, const int atomCharge
 						, const size_t bondValenceSum
 						, const size_t bondNum
 						, const size_t bondNumAromatic
 						, const size_t bondNumProton
 						);


		  /*! Adds explicit "H" atoms according to the valence of atom nodes
		   * and protons within complex atom labels.
		   *
		   * TODO: NOTE, this implementation is still experimental and is not
		   * necessarily correct!
 		   *
		   * @param mol the molecule graph to fill
		   */
		static
		void
		fillProtons( Molecule& mol );


 		  /*!
 		   * Converts a given molecule in Chemical Markup Language (CML) format.
 		   *
 		   * See http://cml.sourceforge.net/ for further details on the format.
 		   *
 		   * @param mol the molecule to represent
 		   * @return the CML string representation of the molecule
 		   */
 		static
 		std::string
 		convertCML( const Molecule& mol );


 		  /*!
 		   * Converts a given molecule in Chemical Markup Language (CML) format.
 		   *
 		   * See http://cml.sourceforge.net/ for further details on the format.
 		   *
 		   * @param mol the molecule to represent
 		   * @param stream the out stream to add the CML representation to
 		   * @return the altered stream filled with the CML string
 		   *         representation of the molecule
 		   */
 		static
 		std::ostream &
 		convertCML( const Molecule& mol, std::ostream & stream );


 		  /*!
 		   * Replaces all nodes with group labels with the according molecule
 		   * component if present in groups. If the group label is unknown, an
 		   * exception is raised.
 		   *
 		   * NOTE: The node replacement is NOT recursive to avoid infinite
 		   * replacement chains.
 		   *
 		   * @param mol the molecule to alter
 		   * @param groups the list of known groups
 		   *
 		   * @thrown std::runtime_error in case a group label is unknown
 		   */
 		static
 		void
 		insertGroups( ggl::chem::Molecule & mol
 					, const GroupMap & groups
 					) throw(std::runtime_error);


 		  /*!
 		   * Iteratively identifies groups within the molecule and replaces them
 		   * by the according group label. The replacement is done in decreasing
 		   * group size, ie. first the largest groups are compressed.
 		   *
 		   * @param mol the molecule to alter
 		   * @param groups the list of known groups to be introduced if found
 		   *
 		   */
 		static
 		void
 		compressGroups( ggl::chem::Molecule & mol
						, const GroupMap & groups
						);

 		  /*!
 		   * Checks whether or not a given atom node label is a group label.
 		   * @param nodeLabel the node label of interest
 		   * @return true if the given node label is a group label; false
 		   *         otherwise
 		   */
 		static
 		bool
 		isGroupLabel( const std::string & nodeLabel );


 		  /*!
 		   * Extracts the group identifier from the (complex) node label.
 		   * @param nodeLabel the node label that contains a group label
 		   * @return the group label or an empty string if no group label is
 		   *         present
 		   */
 		static
 		std::string
 		getGroupLabel( const std::string & nodeLabel );

	protected:
		
		  //! holds all digits for string parsing
		static const std::string DIGIT;
		  //! holds all white space characters for string parsing
		static const std::string WHITESPACE;

		  //! the stored atom label to data mapping
 		static AtomDataMap atomData;

 		  //! maps aromatic labels onto their non-aromatic form and vice versa
 		static AromaticSwapMap aromaticSwapMap;

		  //! the stored bond label to data mapping
 		static BondDataMap bondData;

 		
 		  /*! Checks if all node label of the given molecule are supported by
 		   * the current SMILESwriter.
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff all node label are ok, false otherwise
 		   */
 		static
 		bool
 		checkAtomLabel( const Molecule & mol );
 		
 		  /*! Checks if any node label is equal to AtomLabelWildcard. Since this
 		   * is currently not supported the check returns false if any occurence
 		   * is found.
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff no atom label equals AtomLabelWildcard;
 		   *         false otherwise
 		   */
 		static
 		bool
 		checkAtomLabelWildcard( const Molecule & mol );

 		  /*! Checks if any complex node label contains implicit H atoms.
 		   * This is currently not supported by the library
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff no atom label equals AtomLabelWildcard;
 		   *         false otherwise
 		   */
 		static
 		bool
 		checkAtomComplexWithH( const Molecule & mol );

 		  /*! Checks if any atom shows an inconsistent electron distribution,
 		   * ie. [valence+charge != bondValenceSum+aromAdd] where aromAdd == 1
 		   * if the atom is part of an aromatic ring and 0 otherwise.
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff no atom shows inconsistent electron distribution
 		   *         false otherwise
 		   */
 		static
 		bool
 		checkAtomValence( const Molecule & mol );

 		  /*! Checks if all bond label of the given molecule are supported by
 		   * the current SMILESwriter.
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff all bond labels are ok; false otherwise
 		   */
 		static
 		bool
 		checkBondLabel( const Molecule & mol );

 		  /*! Checks if an atom forms a bond with itself, i.e. a ring bond.
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff all bond labels are ok; false otherwise
 		   */
 		static
 		bool
 		checkBondLoop( const Molecule & mol );

 		  /*! Checks if the molecule shows more than one connected component.
 		   *
 		   * @param mol the molecule to check
 		   * @return true iff only the molecule is connected; false otherwise
 		   */
 		static
 		bool
 		checkNonConnected( const Molecule & mol );



	};
	
	

  } // namespace chem
}  // namespace ggl

#include "ggl/chem/MoleculeUtil.icc"

#endif /*MOLECULEUTIL_HH_*/
