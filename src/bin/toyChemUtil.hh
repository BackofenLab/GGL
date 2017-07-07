#ifndef TOYCHEMUTIL_HH_
#define TOYCHEMUTIL_HH_

		
//////////////////////////////////////////////////////////////////////////

#include <string>

#include <boost/version.hpp>

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

#include <ggl/chem/Molecule.hh>
#include <ggl/chem/MoleculeComponent.hh>
#include <ggl/chem/MoleculeUtil.hh>


  //////////////////////////////////////////////////////////////////////////
 /*! The used container to store SMILES and their corresponding molecule graph
  * objects.
  */
typedef
#if HAVE_UNORDERED_MAP > 0
	std::unordered_map<std::string, ggl::chem::Molecule *>
#elif HAVE_TR1_UNORDERED_MAP > 0
	std::tr1::unordered_map<std::string, ggl::chem::Molecule *>
#elif HAVE_GNU_HASH_MAP > 0
	__gnu_cxx::hash_map< std::string, ggl::chem::Molecule*, sgm::hash_string >
#else
	std::map< std::string, ggl::chem::Molecule* >
#endif
		 SMILES_container;
  //////////////////////////////////////////////////////////////////////////

		
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <ggl/chem/ChemRule.hh>
#include <ggl/chem/ChemRuleGraph.hh>
#include <ggl/chem/AromaticityPerception.hh>
		
  //////////////////////////////////////////////////////////////////////////
  /*! The used container to store the LeftSidePattern of Rules grouped by the
   * number of components of the LeftSidePattern.
   */
typedef
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map< size_t, std::vector< ggl::chem::LeftSidePattern* > >
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map< size_t, std::vector< ggl::chem::LeftSidePattern* > >
#elif HAVE_GNU_HASH_MAP > 0
		__gnu_cxx::hash_map< size_t, std::vector< ggl::chem::LeftSidePattern* > >
#else
		std::map< size_t, std::vector< ggl::chem::LeftSidePattern* > >
#endif
		RulePatternMap;
  //////////////////////////////////////////////////////////////////////////
		 
			
//////////////////////////////////////////////////////////////////////////

		 
#include <exception>
#include <string>
		 
#ifndef ARGEXCEPTION_
#define ARGEXCEPTION_
	
	  /*! Exception class for exeptions thrown during argument and input parsing.
	   */
	class ArgException : public std::exception {
	public:
		  //! the error message
		std::string errorMsg;
		ArgException( std::string errorMsg_ ) : errorMsg(errorMsg_) {}
		virtual ~ArgException() throw (){}
		virtual const char* what() const throw() {
			return errorMsg.c_str();
		}
	};

#endif
	
//////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
	
#include <ggl/RuleGraph.hh>
#include <ggl/chem/Molecule.hh>
#include "ggl/chem/ChemRule.hh"


//////////////////////////////////////////////////////////////////////////
	
	
	/*!
	 * Cuts a ':' separated list of files and stores all non-empty list elements
	 * in the returned vector.
	 * 
	 * @param list the list of file names
	 * @return the list of non-empty list elements
	 */
	std::vector< std::string >
	splitToFileNames( const std::string & list );
	

//////////////////////////////////////////////////////////////////////////

	
	/*!
	 * Parses the given input sources for rules and adds each rule to the given
	 * container.
	 * 
	 * @param inSource either STDIN for getting input from standard input, or a
	 *       ':' separated list of file names to read from
	 * @param toFill the container to add the rules to
	 * @param groups molecule groups IDs and the according molecule components
	 *       to insert instead
	 * @param isCompactGML whether or not the stream to read contains compacted
	 *       rule GML notation or not
	 */
	void
	parseRules(	const std::string & inSource
				, std::vector<ggl::chem::ChemRule> & toFill
				, const ggl::chem::GroupMap & groups
				, const bool compactGML = false) throw(std::exception);

	
//////////////////////////////////////////////////////////////////////////


	/*!
	 * Parses the given input stream for rules and adds each rule to the given
	 * container.
	 * 
	 * Each rule has to be headed by a comment line starting with a '#' 
	 * character that specifies the name of the rule.
	 * 
	 * @param in the input stream to read from
	 * @param toFill the container to add the rules to
	 * @param groups molecule groups IDs and the according molecule components
	 *       to insert instead
	 * @param linesRead the number of lines already read from input (needed for 
	 *                  error reporting)
	 * 
	 */
	void
	parseRules(	std::istream & in
				, std::vector<ggl::chem::ChemRule> & toFill
				, const ggl::chem::GroupMap & groups
				, const size_t linesRead ) throw(std::exception);  

//////////////////////////////////////////////////////////////////////////

	
	/*!
	 * Parses the given input sources for SMILES strings and adds each parsed 
	 * molecule to the given container.
	 * 
	 * Each line in each input source should contain only ONE SMILES string. 
	 * Leading and tailing whitespaces are ignored.
	 * 
	 * @param inSource either STDIN for getting input from standard input, or a
	 *       ':'-separated list of file names to read from
	 * @param toFill the container to add the found SMILES and molecules to
	 * @param groups molecule groups IDs and the according molecule components
	 *       to insert instead
	 * @param setNextAtomClass if > 0, atoms will be numbered starting from the
	 *       given value and the numbering will be added as class label;
	 *       replace class labels will be reported
	 * @return if setNextAtomClass > 0, the successive atom class label that was
	 *       not used yet (for iterative calls); 0 otherwise
	 */
	size_t
	parseSMILES(	const std::string & inSource
					, SMILES_container & toFill
					, const ggl::chem::GroupMap & groups
					, const size_t setNextAtomClass = 0
					) throw(std::exception);

	
//////////////////////////////////////////////////////////////////////////

	
	/*!
	 * Parses SMILES strings from stream. Each line should contain only 
	 * ONE SMILES string. Leading and tailing whitespaces are ignored.
	 * 
	 * @param in the stream to read from
	 * @param toFill the container to add the found SMILES and molecules to
	 * @param linesRead the number of lines already read from input (needed for 
	 *                  error reporting)
	 * @param groups molecule groups IDs and the according molecule components
	 *       to insert instead
	 * @param setNextAtomClass if > 0, atoms will be numbered starting from the
	 *       given value and the numbering will be added as class label;
	 *       replace class labels will be reported
	 * @return if setNextAtomClass > 0, the successive atom class label that was
	 *       not used yet (for iterative calls); 0 otherwise
	 */
	size_t
	parseSMILES(	std::istream & in
					, SMILES_container & toFill
					, const size_t linesRead
					, const ggl::chem::GroupMap & groups
					, const size_t setNextAtomClass = 0
					) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

	
	/*!
	 * Parses the given input source for graphs in GML format and adds each 
	 * graph to the given container.
	 * 
	 * @param inSource either STDIN for getting input from standard input, or a
	 *       ':' separated list of file names to read from
	 * @param toFill the container to add the found SMILES and molecules to
	 * @param groups molecule groups IDs and the according molecule components
	 *       to insert instead
	 * @param setNextAtomClass if > 0, atoms will be numbered starting from the
	 *       given value and the numbering will be added as class label;
	 * @return if setNextAtomClass > 0, the successive atom class label that was
	 *       not used yet (for iterative calls); 0 otherwise
	 */
	size_t
	parseMolGML(	const std::string & inSource
					, SMILES_container & toFill
					, const ggl::chem::GroupMap & groups
					, const size_t setNextAtomClass = 0
					) throw(std::exception);

	
//////////////////////////////////////////////////////////////////////////
	

	/*!
	 * Parses the given input stream for graphs in GML format and adds each 
	 * graph to the given container.
	 * 
	 * Each graph has to be headed by a comment line starting with a '#' 
	 * character that specifies the name of the graph.
	 * 
	 * @param in the input stream to read from
	 * @param toFill the container to add the found SMILES and molecules to
	 * @param linesRead the number of lines already read from input (needed for 
	 *                  error reporting)
	 * @param groups molecule groups IDs and the according molecule components
	 *       to insert instead
	 * @param pruneProtons if true, all protons are removed from the molecule
	 *           graph before generating the SMILES string; otherwise they are
	 *           compressed into complex labels from adjacent non-proton atoms
	 * @param reportConversion if true, the molecule graph resulting from the
	 *           proton pruning/compression is reported to STDOUT; otherwise
	 *           no reporting is done.
	 * @param error if non-NULL, each error is reported to this stream and no
	 *           exception is thrown; otherwise errors are handled reported via
	 *           exceptions.
	 * @param correctNodeByBonds if true, aromaticity of node labels is correct
	 *           according to the presence of adjacent aromatic bonds.
	 * @param setNextAtomClass if > 0, atoms will be numbered starting from the
	 *       given value and the numbering will be added as class label;
	 * @return if setNextAtomClass > 0, the successive atom class label that was
	 *       not used yet (for iterative calls); 0 otherwise
	 */
	size_t
	parseMolGML(	std::istream & in
					, SMILES_container & toFill
					, const size_t linesRead
					, const ggl::chem::GroupMap & groups
					, const size_t setNextAtomClass = 0
					, const bool pruneProtons = false
					, const bool reportConversion = false
					, std::ostream *error = NULL
					, const bool correctNodeByBonds = false
					) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

	 /*!
	 * Parses the given input source for molecule components in GML format and
	 * adds each component to the given container.
	 *
	 * @param inSource either STDIN for getting input from standard input, or a
	 *       ':' separated list of file names to read from
	 * @param toFill the container to add the found molecule components and
	 *       according identifier to
	  */
	void
	parseGroups( const std::string& inSource, ggl::chem::GroupMap& toFill ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

	/*!
	 * Parses GML molecule group encodings from stream.
	 *
	 * @param in the stream to read from
	 * @param toFill the inserter to add the found SMILES to
	 *
	 */
	void
	parseGroups(	std::istream & in
					, ggl::chem::GroupMap & toFill ) throw(std::exception);

//////////////////////////////////////////////////////////////////////////

	void
	printRule( std::ostream& out, ggl::chem::ChemRule & rule );

//////////////////////////////////////////////////////////////////////////

	void
	printRules( std::ostream& out, std::vector<ggl::chem::ChemRule> & rules );

//////////////////////////////////////////////////////////////////////////

	void
	printCopyAndPaste( std::ostream & out, const ggl::chem::ChemRule & rule );

//////////////////////////////////////////////////////////////////////////

	void
	printConstraints( std::ostream & out, const ggl::chem::ChemRule & rule );

//////////////////////////////////////////////////////////////////////////

	void
	printSMILES( std::ostream& out, SMILES_container & smiles );

//////////////////////////////////////////////////////////////////////////
	
	void
	printGraph( std::ostream& out, const sgm::Graph_Interface& g );

//////////////////////////////////////////////////////////////////////////

	void
	printGroups( std::ostream& out, const ggl::chem::GroupMap & groups );

//////////////////////////////////////////////////////////////////////////
	
	struct CompareStringPointer {
		bool operator()(const std::string * const& p1, const std::string * const& p2 ) const {
			return p1->size() < p2->size()
					|| (p1->size() == p2->size() && *p1 < *p2);
		}
	};

	struct CompareString {
		bool operator()(const std::string & p1, const std::string & p2 ) const {
			return p1.size() < p2.size()
					|| (p1.size() == p2.size() && p1 < p2);
		}
	};

//////////////////////////////////////////////////////////////////////////

#include <vector>

#include <sgm/SubGraphMatching.hh>
#include <sgm/Match_Reporter.hh>
#include <sgm/Pattern_Automorphism.hh>
#include <ggl/RuleGraph.hh>
#include <ggl/chem/Molecule.hh>
	
	
//////////////////////////////////////////////////////////////////////////


	/*!
	 * Performs the application of one rule to all given targets by calling the 
	 * recursive handling by singleRuleApplicationRec.
	 * 
	 * @param sgm the SubGraphMatching implementation to use
	 * @param rulePattern the LeftSidePattern of a Rule to be searched and 
	 *                    applied
	 * @param allTargets the list of all Molecules to apply the rule onto
	 * @param mrApplyRule the match reporter that applies the Rule of 
	 *                    rulePattern
	 * @param ruleSymmetry the symmetry definition of the rulePattern to avoid
	 *                     match reporting overhead
	 * @param ignoreAtomClassLabel whether or not the atom class labels are to
	 *        be ignored for the rule pattern matching, i.e. all the atom
	 *        labels from all initialMolecules are reduced to the atom label
	 *        prefix substring up to the class label separator ':' (if present)
	 * @param noRedundantMolecules whether or not the same molecule is allowed
	 *        to be part of the target molecule compilation for multi-component
	 *        left-side rule pattern compilations
	 * @return true if this component and all of the remaining were present to 
	 * allow for an early abortion if one of the component cannot be matched
	 */
	bool 
	singleRuleApplication(	sgm::SubGraphMatching& sgm
							, const ggl::chem::LeftSidePattern& rulePattern
							, const SMILES_container & allTargets
							, sgm::Match_Reporter& mrApplyRule
							, const sgm::Pattern_Automorphism* ruleSymmetry=NULL
							, const bool ignoreAtomClassLabel = false
							, const bool noRedundantMolecules = false
							);

//////////////////////////////////////////////////////////////////////////

	/*!
	 * Performs the application of one rule to all given targets by generating 
	 * all possible combinations of targets for each component to keep the 
	 * search space of the SubGraphMatching implementation as small as possible.
	 * 
	 * @param sgm the SubGraphMatching implementation to use
	 * @param rulePattern the LeftSidePattern of a Rule to be searched and 
	 *                    applied
	 * @param ruleComponent the current component number to check within the 
	 *                      recursion
	 * @param allTargets the list of all Molecules to apply the rule onto
	 * @param curTargets the list of one Molecule per Rule pattern component to
	 *                   be used as the graph to match against
	 * @param mrApplyRule the match reporter that applies the Rule of 
	 *                    rulePattern
	 * @param ruleSymmetry the symmetry definition of the rulePattern to avoid
	 *                     match reporting overhead
	 * @param ignoreAtomClassLabel whether or not the atom class labels are to
	 *        be ignored for the rule pattern matching, i.e. all the atom
	 *        labels from all initialMolecules are reduced to the atom label
	 *        prefix substring up to the class label separator ':' (if present)
	 * @param noRedundantMolecules whether or not the same molecule is allowed
	 *        to be part of the target molecule compilation for multi-component
	 *        left-side rule pattern compilations
	 * @return true if this component and all of the remaining were present to 
	 * allow for an early abortion if one of the component cannot be matched
	 */
	bool 
	singleRuleApplicationRec(	sgm::SubGraphMatching& sgm
						, const ggl::chem::LeftSidePattern& rulePattern
						, const size_t ruleComponent
						, const SMILES_container & allTargets
						, std::vector< const ggl::chem::Molecule* >& curTargets
						, sgm::Match_Reporter& mrApplyRule
						, const sgm::Pattern_Automorphism* ruleSymmetry = NULL
						, const bool ignoreAtomClassLabel = false
						, const bool noRedundantMolecules = false
						);

//////////////////////////////////////////////////////////////////////////

	
//////////////////////////////////////////////////////////////////////////


	/*!
	 * Performs the application of one rule to all given targets by calling the 
	 * recursive handling by singleRuleApplicationRec.
	 * 
	 * @param sgm the SubGraphMatching implementation to use
	 * @param rulePattern the LeftSidePattern of a Rule to be searched and 
	 *                    applied
	 * @param newTargets the list of all Molecules to apply the rule onto
	 * @param oldTargets the list of all Molecules to apply the rule onto but
	 *                    that have been matched by the rule before without the
	 *                    newTargets
	 * @param mrApplyRule the match reporter that applies the Rule of 
	 *                    rulePattern
	 * @param ruleSymmetry the symmetry definition of the rulePattern to avoid
	 *                     match reporting overhead
	 * @param ignoreAtomClassLabel whether or not the atom class labels are to
	 *        be ignored for the rule pattern matching, i.e. all the atom
	 *        labels from all initialMolecules are reduced to the atom label
	 *        prefix substring up to the class label separator ':' (if present)
	 * @param noRedundantMolecules whether or not the same molecule is allowed
	 *        to be part of the target molecule compilation for multi-component
	 *        left-side rule pattern compilations
	 * @return true if this component and all of the remaining were present to 
	 * allow for an early abortion if one of the component cannot be matched
	 */
	bool 
	singleRuleApplication(	sgm::SubGraphMatching& sgm
							, const ggl::chem::LeftSidePattern& rulePattern
							, const SMILES_container & newTargets
							, const SMILES_container & oldTargets
							, sgm::Match_Reporter& mrApplyRule
							, const sgm::Pattern_Automorphism* ruleSymmetry=NULL
							, const bool ignoreAtomClassLabel = false
							, const bool noRedundantMolecules = false
							);

//////////////////////////////////////////////////////////////////////////

	/*!
	 * Performs the application of one rule to all given targets by generating 
	 * all possible combinations of targets for each component to keep the 
	 * search space of the SubGraphMatching implementation as small as possible.
	 * 
	 * @param sgm the SubGraphMatching implementation to use
	 * @param rulePattern the LeftSidePattern of a Rule to be searched and 
	 *                    applied
	 * @param ruleComponent the current component number to check within the 
	 *                      recursion
	 * @param newTargets the list of all Molecules to apply the rule onto
	 * @param curTargets the list of one Molecule per Rule pattern component to
	 *                    be used as the graph to match against
	 * @param atLeastOneNewAdded TRUE if curTargets contains at least one 
	 *                    molecule from the newTargets container
	 * @param mrApplyRule the match reporter that applies the Rule of 
	 *                    rulePattern
	 * @param ruleSymmetry the symmetry definition of the rulePattern to avoid
	 *                     match reporting overhead
	 * @param ignoreAtomClassLabel whether or not the atom class labels are to
	 *        be ignored for the rule pattern matching, i.e. all the atom
	 *        labels from all initialMolecules are reduced to the atom label
	 *        prefix substring up to the class label separator ':' (if present)
	 * @param noRedundantMolecules whether or not the same molecule is allowed
	 *        to be part of the target molecule compilation for multi-component
	 *        left-side rule pattern compilations
	 * @return 0 if this component and all of the remaining were present
	 *         1 if one of the remaining components cannot be matched to 
	 *           allow for an early abortion
	 *         2 if no newTarget element could be matched on the remaining 
	 *           components but none was matched so far 
	 */
	int 
	singleRuleApplicationRec(	sgm::SubGraphMatching& sgm
						, const ggl::chem::LeftSidePattern& rulePattern
						, const size_t ruleComponent
						, const SMILES_container & newTargets
						, const SMILES_container & oldTargets
						, const bool atLeastOneNewAdded
						, std::vector< const ggl::chem::Molecule* >& curTargets
						, sgm::Match_Reporter& mrApplyRule
						, const sgm::Pattern_Automorphism* ruleSymmetry = NULL
						, const bool ignoreAtomClassLabel = false
						, const bool noRedundantMolecules = false
						);

//////////////////////////////////////////////////////////////////////////

#include "ggl/chem/MR_Reactions.hh"

//////////////////////////////////////////////////////////////////////////

	
	 /*! Applies all rules onto a set of initial molecules. The resulting 
	  * molecules from the applications as well as the reaction information is
	  * written to provided containers.
	  * 
	  * @param rules (IN) the left side pattern of the rules to apply 
	  * @param initialMolecules (IN) the molecules the rules are applied on
	  * @param producedMolecules (OUT) the container where the molecules 
	  *         produced by the rule application are added
	  * @param producedReactions (OUT) the container where the reaction 
	  *         information of the rule application is stored in
	  * @param rateCalc the rate calculation object to use or NULL if none to
	  *         apply
	  * @param aromaticity the aromaticity perception instance to be used to
	  *         correct product molecules generated by rule applications
	  * @param ignoreAtomClassLabel whether or not the atom class labels are to
	  *        be ignored for the rule pattern matching, i.e. all the atom
	  *        labels from all initialMolecules are reduced to the atom label
	  *        prefix substring up to the class label separator ':' (if present)
	  * @param noRedundantMolecules whether or not the same molecule is allowed
	  *        to be part of the target molecule compilation for multi-component
	  *        left-side rule pattern compilations
	  */
	void
	applyRules( const RulePatternMap & rules
				, const SMILES_container & initialMolecules
				, SMILES_container & producedMolecules
				, ggl::chem::MR_Reactions::Reaction_Container& producedReactions
				, const ggl::chem::ReactionRateCalculation * rateCalc
				, const ggl::chem::AromaticityPerception & aromaticity
				, const bool ignoreAtomClassLabel = false
				, const bool noRedundantMolecules = false
		);
	
	
//////////////////////////////////////////////////////////////////////////


	
	 /*! Applies all rules onto a set of initial molecules. The resulting 
	  * molecules from the applications as well as the reaction information is
	  * written to provided containers.
	  * 
	  * @param rules (IN) the left side pattern of the rules to apply 
	  * @param oldMolecules (IN) molecules to that the rules have been applied
	  *        already (e.g. in last iteration), such that no rule application 
	  *        to this set of rules only is done
	  * @param newMolecules (IN) molecules the rules were not applied on 
	  *        already such that at least one of these molecules will be within
	  *        a single rule application target
	  * @param producedMolecules (OUT) the container where the molecules 
	  *         produced by the rule application are added
	  * @param producedReactions (OUT) the container where the reaction 
	  *         information of the rule application is stored in
	  * @param rateCalc the rate calculation object to use or NULL if none to
	  *         apply
	  * @param allowAllIntra (IN) whether or not the application of 
	  *         multicomponent rules is allowed on a single molecule only
	  * @param aromaticity the aromaticity perception instance to be used to
	  *         correct product molecules generated by rule applications
	  * @param ignoreAtomClassLabel whether or not the atom class labels are to
	  *        be ignored for the rule pattern matching, i.e. all the atom
	  *        labels from all initialMolecules are reduced to the atom label
	  *        prefix substring up to the class label separator ':' (if present)
	  * @param noRedundantMolecules whether or not the same molecule is allowed
	  *        to be part of the target molecule compilation for multi-component
	  *        left-side rule pattern compilations
	  */
	void
	applyRules( const RulePatternMap & rules
				, const SMILES_container & oldMolecules
				, const SMILES_container & newMolecules
				, SMILES_container & producedMolecules
				, ggl::chem::MR_Reactions::Reaction_Container& producedReactions
				, const ggl::chem::ReactionRateCalculation * rateCalc
				, const bool allowAllIntra
				, const ggl::chem::AromaticityPerception & aromaticity
				, const bool ignoreAtomClassLabel = false
				, const bool noRedundantMolecules = false
			);
	
	
//////////////////////////////////////////////////////////////////////////

	
//	 /*! Applies all rules onto a set of initial molecules. The resulting 
//	  * molecules from the applications are
//	  * written to provided containers.
//	  * 
//	  * @param rules (IN) the left side pattern of the rules to apply 
//	  * @param oldMolecules (IN) molecules to that the rules have been applied
//	  *        already (e.g. in last iteration), such that no rule application 
//	  *        to this set of rules only is done
//	  * @param newMolecules (IN) molecules the rules were not applied on 
//	  *        already such that at least one of these molecules will be within
//	  *        a single rule application target
//	  * @param producedMolecules (OUT) the container where the molecules 
//	  *         produced by the rule application are added
//	  */
//	void
//	applyRules( const RulePatternMap & rules
//				, const SMILES_container & oldMolecules
//				, const SMILES_container & newMolecules
//				, SMILES_container & producedMolecules
//			);
//	
	
//////////////////////////////////////////////////////////////////////////

	/**
	 * adds an atom class information to all atoms without using consecutive
	 * numbers starting at nextClassLabelOption
	 *
	 * @param mol the molecule to correct
	 * @param setNextAtomClass atoms without class label will be numbered
	 *       starting from the given value and the numbering will be added
	 *       as class label;
	 * @param overwriteExistingClass whether or not to overwrite existing class
	 *       labels; if true, a warning is produced with the new class label
	 *       mapping
	 * @return the successive atom class label that was
	 *       not used yet (for iterative calls)
	 */
	size_t
	setAtomClass( ggl::chem::Molecule & mol
					, const size_t nextClassLabelOption
					, const bool overwriteExistingClass );

//////////////////////////////////////////////////////////////////////////


	/*!
	 *
	 * @param setNextAtomClass if > 0, atoms will be numbered starting from the
	 *       given value and the numbering will be added as class label;
	 * @return if setNextAtomClass > 0, the successive atom class label that was
	 *       not used yet (for iterative calls); 0 otherwise
	 */
	size_t
	correctInputMolecules( const SMILES_container & targetSmiles
							, SMILES_container & producedSmiles
							, ggl::chem::AromaticityPerception * aromaticity_full
							, const bool protonFilling = true
							, const size_t setNextAtomClass = 0
							);

//////////////////////////////////////////////////////////////////////////



#endif /*TOYCHEMUTIL_HH_*/
