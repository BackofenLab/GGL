#ifndef GGL_CHEM_MR_REACTIONS_HH_
#define GGL_CHEM_MR_REACTIONS_HH_

#include <vector>
#include <set>
#include <cassert>

#include "sgm/Match_Reporter.hh"
#include "sgm/Graph_boostV_p.hh"

#include "ggl/RuleGraph.hh"
#include "ggl/MR_ApplyRule.hh"
#include "ggl/MR_ApplyRule_NodeLabelPrefix.hh"

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/SMILESwriter.hh"
#include "ggl/chem/GS_MolCheck.hh"
#include "ggl/chem/GS_SMILES.hh"
#include "ggl/chem/Reaction.hh"
#include "ggl/chem/ReactionTransitionState.hh"
#include "ggl/chem/ReactionRateCalculation.hh"
#include "ggl/chem/AP_NSPDK.hh"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include "sgm/HashMap.hh"


#if HAVE_UNORDERED_MAP > 0
	#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP > 0
	#include <tr1/unordered_map>
#elif HAVE_GNU_HASH_MAP > 0
	#include <ext/hash_map>
#else
	#include <map>
#endif
	
	

namespace ggl {
 namespace chem {

  //////////////////////////////////////////////////////////////////////////////
  //! Type for mapping of SMILES strings to their graph 
  //! implementation
typedef 
#if HAVE_UNORDERED_MAP > 0
	std::unordered_map< std::string, ggl::chem::Molecule*>
#elif HAVE_TR1_UNORDERED_MAP > 0
	std::tr1::unordered_map< std::string, ggl::chem::Molecule*>
#elif HAVE_GNU_HASH_MAP > 0
	__gnu_cxx::hash_map< std::string, ggl::chem::Molecule*, sgm::hash_string >
#else
	std::map< std::string, ggl::chem::Molecule* >
#endif
	Smiles2GraphMap;
  //////////////////////////////////////////////////////////////////////////////

	
  //////////////////////////////////////////////////////////////////////////////
  //! Type for mapping of graph implementations to SMILES strings
typedef	
#if HAVE_UNORDERED_MAP > 0
	std::unordered_map< const ggl::chem::Molecule*, std::string >
#elif HAVE_TR1_UNORDERED_MAP > 0
	std::tr1::unordered_map< const ggl::chem::Molecule*, std::string >
#elif HAVE_GNU_HASH_MAP > 0
	__gnu_cxx::hash_map< const ggl::chem::Molecule*, std::string >
#else
	std::map< const ggl::chem::Molecule*, std::string >
#endif
	Graph2SmilesMap;
  //////////////////////////////////////////////////////////////////////////////
	
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


	  /*! @brief Reaction match reporter
	   *
	   * A Match_Reporter that generates a Reaction object from a rule
	   * application.
	   * 
	   * NOTE : this class currently only works for target graphs that are an
	   * object of Molecule_GraphV !!! Furthermore, it assumes that the
	   * elements (pointer) from the covered target graph are from the here
	   * maintained storage, to enable a fast lookup of their SMILES !!!
	   * Finaly, SMILES are generated using an instance of SMILESwriter !!!
	   * 
	   * NOTE : objects of this class are not safe to copy due to Ruler member!
	   *
	   *  @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */ 
	class MR_Reactions
	  : public sgm::Match_Reporter
	  	, public GS_SMILES_MOLp< Smiles2GraphMap >
	{
		
	public:
		
		  //! a container that stores Reaction objects
		typedef std::set< Reaction > Reaction_Container;
		
	protected:
		
		  //! the SMILESwriter to use for target-to-smiles conversion
		typedef SMILESwriter Smiler;
										
		  //! access to the superclass of this class
		typedef GS_SMILES_MOLp< Smiles2GraphMap > SuperClass;
		
	protected:
		
		  //! The storage of known SMILES to the graph implementations.
		using SuperClass::smiles2molCheckOnly1;
		  //! An alternative second storage of known SMILES to the graph 
		  //! implementations. 
		using SuperClass::smiles2molCheckOnly2;
		  //! The storage of new SMILES to the graph implementations. 
		  //! The pointer are 'new' allocated and have to be destroyed outside!
		using SuperClass::smiles2mol;
		  //! A vice versa lookup table for smiles2graph
		Graph2SmilesMap graph2smiles;
		
		  //! The container of Reaction objects to fill and extend
		Reaction_Container& reactions;
		
		  //! the current Reaction to fill and to store later
		Reaction curReaction;

		  //! if set to true, than all components of the Rules LeftSidePattern 
		  //! are matched to an own target graph copy if two components map to 
		  //! the same target graph
		const bool addEachComponent;
		
		  //! the checker used to correct all molecules produced by the Ruler;
		  //! it will apply an AP_NSPDK aromaticity perception instance.
		GS_MolCheck molchecker;

		  //! the Match_Reporter utilized to apply the rule and to produce the 
		  //! resulting molecules that will be stored in the reaction
		ggl::MR_ApplyRule * Ruler;
							
		  //! if != NULL, for each new reaction a rate is calculated
		const ReactionRateCalculation* rateCalculation;
		
	public:
		
		  /*! Construction
		   * 
		   * @param smiles2graph the initial container that contains all
		   *        molecules known so far and is fixed
		   * @param newSmiles2graph the container where all new molecules will
		   *        be stored in
		   * @param reactions the Reaction container to fill
		   * @param addEachComponent if set to true, than all components of the
		   *        Rules LeftSidePattern are matched to an own target graph
		   *        copy if two components map to the same target graph
		   * @param rateCalculation The object to use to calculate a reaction
		   *        rate for each new reaction created. If not present (==NULL)
		   *        no reaction rate will be calculated.
		   * @param aromaticity the aromaticity perception class to be used to
		   *        correct rings in the product molecules after the rule
		   *        application was done.
		   * @param keepAtomClass whether or not atom class label information is
		   *        to be preserved during rule rewrite
		   */
		MR_Reactions(	const Smiles2GraphMap& smiles2graph
						, Smiles2GraphMap& newSmiles2graph
						, Reaction_Container& reactions
						, const bool addEachComponent = false
						, const ReactionRateCalculation * rateCalculation = NULL
						, const AromaticityPerception & aromaticity =
								AP_NSPDK("Marvin:general:2013")
						, const bool keepAtomClass = false );
		

		  /*! Construction
		   * 
		   * @param smiles2graph the initial container that contains all
		   *        molecules known so far and which is not changed
		   * @param smiles2graph2 second container of the initial
		   *        molecules known so far 
		   * @param newSmiles2graph the container where all new molecules will
		   *        be stored in
		   * @param reactions the Reaction container to fill
		   * @param addEachComponent if set to true, than all components of the
		   *        Rules LeftSidePattern are matched to an own target graph
		   *        copy if two components map to the same target graph,
		   * @param rateCalculation The object to use to calculate a reaction
		   *        rate for each new reaction created. If not present (==NULL)
		   *        no reaction rate will be calculated.
		   * @param aromaticity the aromaticity perception class to be used to
		   *        correct rings in the product molecules after the rule
		   *        application was done.
		   * @param keepAtomClass whether or not atom class label information is
		   *        to be preserved during rule rewrite
		   */
		MR_Reactions(	const Smiles2GraphMap& smiles2graph
						, const Smiles2GraphMap& smiles2graph2
						, Smiles2GraphMap& newSmiles2graph
						, Reaction_Container& reactions
						, const bool addEachComponent = false
						, const ReactionRateCalculation * rateCalculation = NULL
						, const AromaticityPerception & aromaticity =
								AP_NSPDK("Marvin:general:2013")
						, const bool keepAtomClass = false
						);
		
		  /*! Destruction
		   */
		virtual
		~MR_Reactions();
		
		
		/////////////  sgm::Match_Reporter interface  //////////////////////////
		
		  //! Reports a match. The match is encoded using a vector. The length 
		  //! of the vector corresponds to the number of vertices in the pattern
		  //! and position i encodes the matched position of pattern node i in
		  //! the target graph (an instance of sgm::Graph_boostV_p).
		  //! @param pattern the pattern graph that was searched for. NOTE: HAS
		  //!         TO BE AN INSTANCE OF ggl::chem::LeftSidePattern !!!
		  //! @param target the graph the pattern was found within. NOTE: HAS TO
		  //!         BE AN INSTANCE OF sgm::Graph_boostV_p !!!
		  //! @param match contains the indices of the matched pattern nodes in
		  //! the target graph. match[i] corresponds to the mapping of the ith
		  //! vertex in the pattern graph.
		virtual
		void
		reportHit (	const sgm::Pattern_Interface & pattern,
					const sgm::Graph_Interface & target,
					const sgm::Match & match );
		
		
	protected:
		
		/////////////  GS_SMILES_MOLp interface  //////////////////
		
		/*! Writes a given SMILES and the corresponding molecule graph to the
		 * global storage. Furthermore, the given SMILES is added to the
		 * products of the currently processed ggl::chem::Reaction.
		 *
		 * @param SMILES the SMILES string produced
		 * @param molGraph the molecule graph that corresponds to the given SMILES
		 * @return true it the SMILES was not seen so far; false otherwise
		 */
		virtual
		bool
		insert2map(const std::string & SMILES, const Molecule& molGraph );
				
	};  // MR_Reactions
	
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 } // namespace chem
} // namespace ggl

 // include implementation
#include "ggl/chem/MR_Reactions.icc"

#endif /*MR_REACTIONS_HH_*/
