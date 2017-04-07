#ifndef GGL_CHEM_CHEMRULE_HH_
#define GGL_CHEM_CHEMRULE_HH_

#include "ggl/Rule.hh"
#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeUtil.hh"
#include <map>
#include <iostream>

namespace ggl {
namespace chem {

////////////////////////////////////////////////////////////////////////////////

	  /*! @brief Consistency codes and constants
	   *
	   * A helper class that holds the standard ChemRule constants that are
	   * independent of any template parameter.
	   * 
	   * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class ChemRule_Constants : public Rule_ConCodes {
	public:
		  //! consistency code : everything fine
		using Rule_ConCodes::C_Consistent;
		  //! consistency code : no rule ID is present
		using Rule_ConCodes::C_NoID;

		  //! consistency code : at least one node is created or destroyed
		static const size_t C_NodeInsertDeletion;
		  //! consistency code : at least for one bond too many edges are given
		static const size_t C_EdgeMultiplicity;
		  //! consistency code : an atom label is not SMILES conform
		static const size_t C_AtomLabelInvalid;
		  //! consistency code : a bond label is not SMILES conform
		static const size_t C_BondLabelInvalid;
		  //! consistency code : the valence/electron change of an atom is not
		  //!                    balanced
		static const size_t C_UnbalancedElectrons;
		  //! consistency code : at least one atom label is complex and contains
		  //! implicit H atoms, this is currently not supported
		static const size_t C_AtomComplexWithH;
		  //! consistency code : at least one atom label is complex and contains
		  //! implicit H atoms, this is currently not supported
		static const size_t C_BondLoop;


	};


////////////////////////////////////////////////////////////////////////////////


	  /*! @brief Chemical graph grammar rule
	   *
	   * A description of a graph grammar rule with left and righ rule side and
	   * context description. It describes a chemical reaction as well as its
	   * transition state etc.
	   *
	   * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	   *
	   */
	class ChemRule
	  : public ggl::Rule
	   , public ChemRule_Constants
	{
		
	public:
		/////////////////  PUBLIC  ////////////////////////////////
		
		  //! Access to the super class type
		typedef ggl::Rule SuperClass;
		
		  //! Access to the CoreGraph template argument.
		typedef SuperClass::CoreGraph CoreGraph;
		  //! Access to the NODE_INDEX_PROPERTY template argument.
		typedef SuperClass::NodeIndexProperty NodeIndexProperty;
		  //! Access to the NODE_RULE_CONTEXT_PROPERTY template argument.
		typedef SuperClass::NodeContextProperty NodeContextProperty;
		  //! Access to the EDGE_RULE_CONTEXT_PROPERTY template argument.
		typedef SuperClass::EdgeContextProperty EdgeContextProperty;
		  //! Access to the NODE_LABEL_PROPERTY template argument.
		typedef SuperClass::NodeLabelProperty NodeLabelProperty;
		  //! Access to the NODE_RIGHT_LABEL_PROPERTY template argument.
		typedef SuperClass::NodeRightLabelProperty NodeRightLabelProperty;
		  //! Access to the EDGE_LABEL_PROPERTY template argument.
		typedef SuperClass::EdgeLabelProperty EdgeLabelProperty;
		  //! A container that stores the information of a context side of 
		  //! a graph grammar rule like left, right, and context.
		typedef SuperClass::RuleSide RuleSide;
		  //! Cut-and-paste operation definition
		typedef SuperClass::RuleCnP RuleCnP;
		  //! Cut-and-paste operation container
		typedef SuperClass::CopyAndPasteOperations CopyAndPasteOperations;

		
		/////////////////  CONSISTENCY TYPEDEFS  /////////////////////////
		
		  //! consistency code : everything fine 
		using ChemRule_Constants::C_Consistent;
		  //! consistency code : at least one node is created or destroyed
		using ChemRule_Constants::C_NodeInsertDeletion;
		  //! consistency code : at least for one bond too many edges are given
		using ChemRule_Constants::C_EdgeMultiplicity;
		  //! consistency code : no rule/reaction ID is present
		using ChemRule_Constants::C_NoID;
		  //! consistency code : an atom label is not SMILES conform
		using ChemRule_Constants::C_AtomLabelInvalid;
		  //! consistency code : a bond label is not SMILES conform
		using ChemRule_Constants::C_BondLabelInvalid;
		  //! consistency code : the valence/electron change of an atom is not
		  //!                    balanced
		using ChemRule_Constants::C_UnbalancedElectrons;
		  //! consistency code : an atom forms a bond with itself
		using ChemRule_Constants::C_BondLoop;
		
		/////////////////  CONSTRUCTION  ////////////////////////////////
		
		  /*! Construction of a rule based on the given boost graph based core
		   * description.
		   * @param core the boost graph that describes the graph grammar rule
		   * @param id the identifier or description of the graph grammar rule
		   */
		ChemRule (	const CoreGraph & core
				, const std::string& id = std::string("") );

		  /*! Construction of a rule based on the given boost graph based core
		   * description.
		   * @param core the boost graph that describes the graph grammar rule
		   * @param id the identifier or description of the graph grammar rule
		   * @param constraints the additional constraints for the
		   *        matching of the left side of the rule
		   * @param copyAndPaste the Cut-and-Paste operations to perform
		   */
		ChemRule (	const CoreGraph & core
				, const std::string& id
				, const std::vector< sgm::Pattern_Interface::Match_Constraint* > & constraints
				, const CopyAndPasteOperations & copyAndPaste );

		  /*! Construction of a rule based on the given boost graph based core
		   * description.
		   * @param core the boost graph that describes the graph grammar rule
		   * @param id the identifier or description of the graph grammar rule
		   * @param constraints the additional constraints for the
		   *        matching of the left side of the rule
		   */
		ChemRule (	const CoreGraph & core
				, const std::string& id
				, const std::vector< sgm::Pattern_Interface::Match_Constraint* > & constraints );
		
		  //! Copy construction.
		  //! @param toCopy the ChemRule to make this a copy of
		ChemRule ( const ChemRule & toCopy );
		
		  //! Construction from a normal graph grammar rule.
		  //! @param toCopy the ChemRule to make this a copy of
		ChemRule ( const ggl::Rule & toCopy );

		  //! Destruction
		virtual 
		~ChemRule();
		
		  //! Assignment operator
		  //! @param toCopy the ChemRule to make this a copy of
		  //! @return *this
		ChemRule &
		operator =( const ChemRule & toCopy );
		
		  /*! Checks whether or not this rule is a consistent one. If not, a
		   * combination of according error codes is given.
		   * @return the value is either C_Consistent, or a product of all other
		   *         consistency error codes encountered
		   */
		virtual
		size_t
		isConsistent( void ) const;
		
		  /*!
		   * Writes a description of the consistency status or errors, encoded
		   * in a consistency code produced by isConsistent*(...), to a given
		   * outstream. The function returns whether or not an error occured.
		   *
		   * @param consistencyCode the error code to parse, produced by
		   *           a call to isConsistent(...)
		   * @param errorStream the output stream to write the error decription
		   *           to
		   * @param completeCheck if true: tries to decode the whole
		   *           consistencyCode and reports an error message if this is
		   *           not possible. if false: decodes and reports only the
		   *           known error codes from consistencyCode
		   * @return true if no error is encoded; false otherwise
		   */
		virtual
		bool
		decodeConsistencyStatus(	const size_t consistencyCode
									, std::ostream& errorStream
									, const bool completeCheck = true ) ;


		  /*!
		   * Replaces all nodes with group labels with the according molecule
		   * component if present in groups. If the group label is unknown, the
		   * method returns false; true otherwise.
		   *
		   * NOTE: The node replacement is NOT recursive to avoid infinite
		   * replacement chains.
		   *
		   * @param core the rule core graph to alter
		   * @param groups the list of known groups to be inserted
		   *
		   * @return true if all group labels have been replaced; false
		   *         otherwise
		   *
		   * @throw std::runtime_error if a group is
		   *        to be insert for a non-context node
		   */
		static
		bool
		insertGroups( CoreGraph & core
					, const GroupMap & groups
					) throw(std::runtime_error);



		  /*!
		   * Access to the chemical wildcard string to use for left side
		   * pattern matching
		   * @return the wildcard label MoleculeUtil::AtomLabelWildcard
		   */
		virtual
		const std::string*
		getUsedWildcard(void) const;


		  /*!
		   * Overwrites and disables the wildcard setup. This function does
		   * NOT change the wildcard, since the wildcard is fixed to
		   * MoleculeUtil::AtomLabelWildcard for chemical rules. Thus, this
		   * method is a dummy overwrite.
		   * @param wildcardToUse ignored parameter
		   */
		virtual
		void
		setUsedWildcard( const std::string* const wildcardToUse );



	public:
		/////////////  TRANSITION STATE HANDLING  ///////////////////////////
		
		 /*! Encodes the minimal transition state of a chemical reaction 
		  * that is encoded by a (chemical) graph grammar rule.
		  *
		  * Therein, valence changes (!=0) for each atom are identified and
		  * appended to the class information. All such values show a leading
		  * "000" to make them distinguishable from previously present class
		  * information (which will be the leading number in front of the "000"
		  * separator.
		  *
		  */
		class TransitionState
		{
			
		protected:
			
			  /*! Default construction which is only accessible for subclasses
			   */
			TransitionState();
			
		public:
			/////////////////  PUBLIC  ////////////////////////////////
			
			 /*! Derives a transition state graph from a given graph grammar rule.
			  * 
			  * @param rule the graph grammar rule that encodes the transition state
			  */
			TransitionState( const ChemRule& rule );
			
			 /*! Copy construction
			  * 
			  * @param toCopy the object to make this a copy of
			  */
			TransitionState( const TransitionState& toCopy );
			
			  //! destruction
			virtual 
			~TransitionState();
			
			
			/////////////////  ACCESS MEMBERS ////////////////////////////////
			
			  /*! Allows access to the graph encoding of the transition state.
			   * @return the graph encoding of the transition state
			   */
			const Molecule&
			getGraph( void ) const;
			
			  /*! Allows access to the SMILES encoding of the transition state.
			   * @return the SMILES encoding of the transition state
			   */
			const std::string &
			getSMILES( void ) const;
			
		protected:
			/////////////////  PROTECTED  ////////////////////////////////
			
			  //! the graph encoding of the transition state
			Molecule tState;
			
			  //! SMILES representation of the transition state
			mutable std::string tStateSMILES;
			
		protected:
			
			  //! mapping of vertices to valence change internally used to sum over
			  //! edge valence increases
			typedef std::map< CoreGraph::vertex_descriptor, size_t >
				ValChangeMap;

			 /*! Initializes the transition state graph from a given graph 
			  * grammar rule.
			  * 
			  * @param rule the graph grammar rule that encodes the transition state
			  */
			void
			initializeTransitionState( const CoreGraph& rule );
			
			  /*! Updates a given ValChangeMap : the according entry is increased
			   * with the given valence change value.
			   * 
			   * @param node the entry key to update
			   * @param valChange the value to add to the current value
			   * @param toUpdate the mapping to update within 
			   */
			static
			void
			updateMap(	const CoreGraph::vertex_descriptor& node
						, const size_t valChange
						, ValChangeMap & toUpdate );
			
			  /*! Copies a (modified) rule graph into a molecule object. Only the
			   * primary node and edge label is copied, all other rule information
			   * is ignored.
			   * 
			   * @param rule the rule to copy from
			   * @param mol the molecule to fill / copy to (NOTE: is assumed to be 
			   *            a new/empty graph)
			   */
			static
			void
			copyRuleToMol( const CoreGraph & rule
								, Molecule & mol );
			
		};
		

		  /*! Access to the minimal transition state of the chemical reaction 
		   * covered by this chemical graph grammar rule.
		   * @return the minimal transition state for this rule
		   */ 
		const TransitionState&
		getTransitionState() const;
		
		
	protected:
		/////////////////  PROTECTED  ////////////////////////////////
	
		  //! the boost graph that is encoding the graph grammar rule
		using SuperClass::core;
		  //! the information describing the left side pattern of the Rule
		using SuperClass::LeftSide;
		  //! the information describing the right side (result) of the Rule
		using SuperClass::RightSide;
		  //! the information describing the invariant context of the Rule
		using SuperClass::Context;
		  //! the encoding of the minimal transition state of this reaction
		  //! which is generated on demand (initially not calculated)
		mutable TransitionState* transitionState;
		
	protected:
		
		  /*! Adds no edge constraints for all edges that are right side context
		   * only, to avoid the appearance of parallel edges within molecule
		   * graphs.
		   *
		   * @param toUpdate the ChemRule instance to update
		   */
		static
		void
		addImplicitNoEdgeConstraints( ChemRule & toUpdate );

		  /*! Checks that no node in the graph is created or destroyed, i.e.
		   * if all nodes in the given rule core graph are part of the
		   * context or undergo a label change
		   * 
		   * @param coreGraph the rule core graph to check
		   * @return true if no node is created/destroyed, false otherwise
		   */
		static
		bool
		checkNodeInDel( const CoreGraph & coreGraph );
		
		  /*! Checks that all bonds in the molecule graph are defined by the
		   * correct number of edges. That is either a single context, left, or
		   * right edge, or the combination of a left and right edge. Otherwise,
		   * no single bond between two atoms are defined by the rule.
		   *
		   * @param coreGraph the rule core graph to check
		   * @return true if all bond definitions are correct, false otherwise
		   */
		static
		bool
		checkEdgeInDel( const CoreGraph & coreGraph );

		  /*! Checks if the atom (node) labels are SMILES conform. This is 
		   * checked for input and output of the rule (left/right side).
		   * 
		   * @param coreGraph the rule core graph to check
		   * @return true if all atom/node label are ok, false otherwise
		   */
		static
		bool
		checkAtomLabel( const CoreGraph & coreGraph );
		
		  /*! Checks if any complex node label (in any context) contains
		   * implicit H atoms.
		   * This is currently not supported by the library.
		   *
		   * @param coreGraph the rule core graph to check
		   * @return true if all atom/node labels are ok, false otherwise
		   */
		static
		bool
		checkAtomComplexWithH( const CoreGraph & coreGraph );

		  /*! Checks if the bond (edge) labels are SMILES conform. This is 
		   * checked for input and output of the rule (left/right side).
		   * 
		   * @param coreGraph the rule core graph to check
		   * @return true if all bond/edge label are ok, false otherwise
		   */
		static
		bool
		checkBondLabel( const CoreGraph & coreGraph );
		
		  /*! Checks if an atom forms a bond with itself, i.e. a loop bond.
		   *
		   * @param coreGraph the rule core graph to check
		   * @return true iff all bond labels are ok; false otherwise
		   */
		static
		bool
		checkBondLoop( const CoreGraph & coreGraph );

		  /*! Checks for each atom, if the electron changes of the adjacent 
		   * bonds and the atom label are balanced.
		   * 
		   * NOTE: for aromatic atoms, the an unbalanced electron change of 1
		   * is tolerated due to the flexibility of the involved electrons
		   * within the aromatic ring.
		   *
		   * @param coreGraph the rule core graph to check
		   * @return true if all electron changes are ok, false otherwise
		   */
		static
		bool
		checkElectronChange( const CoreGraph & coreGraph );
		
	};
	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 } // namespace chem
} // namespace ggl

	

// implementation
#include "ggl/chem/ChemRule.icc"

#endif /*GGL_RULE_HH_*/
