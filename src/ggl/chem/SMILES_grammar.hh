
#ifndef GGL_CHEM_SMILES_GRAMMAR_HH_
#define GGL_CHEM_SMILES_GRAMMAR_HH_

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

#include <utility>
#include <vector>
#include <locale>


  // set spirit closure limit if neccessary
#if !defined(BOOST_SPIRIT_CLOSURE_LIMIT)
#define BOOST_SPIRIT_CLOSURE_LIMIT 5
#elif BOOST_SPIRIT_CLOSURE_LIMIT < 5
#error "GGL_CHEM_SMILES_GRAMMAR : BOOST_SPIRIT_CLOSURE_LIMIT too low, has to be at least 5"
#endif

  // set phoenix limit if neccessary
#if !defined(PHOENIX_LIMIT)
#define PHOENIX_LIMIT 5
#elif PHOENIX_LIMIT < 5
#error "GGL_CHEM_SMILES_GRAMMAR : PHOENIX_LIMIT too low, has to be at least 5"
#endif

#include <boost/version.hpp>
#if BOOST_VERSION >= 103800
#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/phoenix1.hpp>
#define NS_BOOSTSPIRIT boost::spirit::classic
#else
#include <boost/spirit.hpp>
#include <boost/spirit/phoenix.hpp>
#define NS_BOOSTSPIRIT boost::spirit
#endif

#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeUtil.hh"

namespace ggl {
  namespace chem {


	  /*! @brief SMILES molecule parser
	   *
	   *  This class defines the rules of the Daylight's (tm) SMILES
	   *  BNF grammar. It allows for the parsing of a SMILES string to generate
	   *  a molecule graph of the encodes molecule. The graph is represented as
	   *  a boost graph (Molecule) and the atom and bond labels will be
	   *  stored in the property_maps of the given PropNodeLabel and
	   *  PropEdgeLabel.
	   *
	   *  BNF grammar of Daylight's SMILES
	   *
	   *  smiles         ::= chain (chain | branch)*
	   *  chain          ::= bond? (simple_atom | complex_atom) ringclosure*
	   *  branch         ::= '(' chain (chain | branch)* ')'
	   *  ringclosure    ::= digit | ('%' digit digit)
	   *  simple_atom    ::= simple_symbol
	   *  complex_atom   ::= '[' isotope? (simple_symbol | complex_symbol | group_symbol) chirality? hcount? charge? name? ']'
	   *  isotope        ::= integer
	   *  simple_symbol  ::= 'Br' | 'Cl' | 'B' | 'c' | etc.
	   *  complex_symbol ::= 's' | 'p' | 'o' | 'Zn' | etc.
	   *  group_symbol   ::= '{' anyChar+ '}'
	   *  chirality      ::= '@' '@'?
	   *  hcount         ::= 'H' integer?
	   *  charge         ::= '+' ('+'* | integer)
	   *                   | '-' ('-'* | integer)
	   *  name           ::= ':' integer
	   *  bond           ::= bond_symbol
	   *  bond_symbol    ::= '-' | '=' | '#' | ':' | etc.
	   *  integer        ::= digit+
	   *  digit          ::= [1-9]
	   *
	   *  NOTE : chirality, and isotope information is currently ignored !
	   *
 	   *  NOTE : Supported atom labels are defined by MoleculeUtil::getAtomData().
 	   *
 	   *  NOTE : Supported bond labels are defined by MoleculeUtil::getBondData().
 	   *
	   *  NOTE further : we allow for an extension of the SMILES encoding:
	   *  complex atoms are allowed to hold a group_symbol ID strings
	   *  enclosed in brackets of the form '{SOMEID}'. They are replaced by
	   *  according group subgraphs if found in the provided group map.
	   *  Otherwise the parsing is aborted.
	   *
	   *  @author Christoph Flamm (c) 2008 http://www.tbi.univie.ac.at/~xtof/
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class SMILES_grammar
	  : public NS_BOOSTSPIRIT::grammar< SMILES_grammar >
	{
	protected:
		
		  //! The boost graph object that is filled to represent the next parsed
		  //! SMILES string.
	    Molecule& g2fill;
	    
	      //! Container that holds group IDs where each matching node has to be
	      //! replaced by the according mapped subgraph
	    const GroupMap * groups;

		  //! Access to the node label property_map of g2fill to set atom labels
	    mutable 
	      boost::property_map<Molecule, PropNodeLabel>::type
	    	nodeLabel;
	    
		  //! Access to the edge label property_map of g2fill to set bond labels
	    mutable 
	      boost::property_map<Molecule, PropEdgeLabel>::type
	    	edgeLabel;
	    
	      //! Adds an atom to the internal molecule graph to fill
	      //! @param label the atom label to set
	    void 
	    addAtom( const std::string& label ) const;
	    
	      //! Adds a bond to the internal molecule graph to fill
	      //! @param atom1 the first bond partner 
	      //! @param atom2 the second bond partner 
	      //! @param label the bond label to set 
	    void
	    addBond( const int atom1, const int atom2, const std::string& label ) const;
	    
	public:
	    
		  //! Constructs the definitions of a Daylight's SMILES grammar to parse
		  //! a SMILES string and to fill the encoded molecule into a given
		  //! boost graph object.
		  //! @param toFill the boost graph object to add nodes and edges to
	    explicit SMILES_grammar( Molecule& toFill );

		  //! Constructs the definitions of a Daylight's SMILES grammar to parse
		  //! a SMILES string and to fill the encoded molecule into a given
		  //! boost graph object.
		  //! @param toFill the boost graph object to add nodes and edges to
	      //! @param groups a container that holds group IDs where
	      //!       each matching node has to be replaced by the according
	      //!       mapped subgraph
	    explicit SMILES_grammar( Molecule& toFill
	    						, const GroupMap & groups );
	    
	      //! Parses a SMILES string and generates a graph representation of the
	      //! molecule
	      //! @param SMILES_string the string to parse
	      //! @return pair.first = the graph encoding of the molecule
	      //!         pair.second = -1 if parsing was successful,
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
		  //! @throw std::invalid_argument in case a check fails
	    static
	    std::pair< Molecule, int >
	    parseSMILES( const std::string & SMILES_string )
	    	throw (std::invalid_argument);
	    
	      //! Parses a SMILES string and generates a graph representation of the
	      //! molecule
	      //! @param SMILES_string the string to parse
	      //! @param groups a container that holds group IDs where
	      //!       each matching node has to be replaced by the according
	      //!       mapped subgraph
	      //! @return pair.first = the graph encoding of the molecule
	      //!         pair.second = -1 if parsing was successful,
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
		  //! @throw std::invalid_argument in case a check fails
	    static
	    std::pair< Molecule, int >
	    parseSMILES( const std::string & SMILES_string, const GroupMap & groups )
	    	throw (std::invalid_argument);

    
		  //! The definition of the SMILES grammar.
	    template <typename ScannerT>
	    struct definition
	    {
	    public:
	    	
	    	  //! Construction of the SMILES BNF grammar rules
	    	definition( SMILES_grammar const& self);
	    	
	    	  //! start parsing
	    	NS_BOOSTSPIRIT::rule<ScannerT> const&
	    	start() const;
	    	
	    protected:

	        char atom2_tmp;
	    	
	         /*!
	          * Dedicated parser for atom labels that can come as simple
	          * symbols and comprise only one characters.
	          */
	        struct simpleSymbol_parser : public NS_BOOSTSPIRIT::char_parser<simpleSymbol_parser>
	        {
	            typedef simpleSymbol_parser self_t;

	            const std::string simpleSymbols;

	            //! construction
	            simpleSymbol_parser()
	             : simpleSymbols("BCNOPSFIHsponcb")
	            {}

	              //! tests whether or not the parsed character is a
	              //! valid and supported atom label
	              //! @param ch the parsed character to test
	              //! @return whether or not @a ch is a valid and supported
	              //! atom label
	            template <typename CharT>
	            bool test(CharT ch) const
	            {
	                return simpleSymbols.find(ch)
	                		!= std::string::npos;
	            }
	        };

	        simpleSymbol_parser const simpleSymbol_p;

	        /*!
	          * Dedicated parser for atom labels comprising only one characters.
	          */
	        struct atom1_parser : public NS_BOOSTSPIRIT::char_parser<atom1_parser>
	        {
	            typedef atom1_parser self_t;

	            //! construction
	            atom1_parser() {}

	              //! tests whether or not the parsed character is a
	              //! valid and supported atom label
	              //! @param ch the parsed character to test
	              //! @return whether or not @a ch is a valid and supported
	              //! atom label
	            template <typename CharT>
	            bool test(CharT ch) const
	            {
	                return MoleculeUtil::getAtomData().find(std::string(1,ch))
	                		!= MoleculeUtil::getAtomData().end();
	            }
	        };

	        atom1_parser const atom1_p;

	          //! first character of a two letter atom label
	        char atom2_firstChar;

	        /*!
	          * Dedicated parser for atom labels comprising only one characters.
	          */
	        struct atom2_parser : public NS_BOOSTSPIRIT::char_parser<atom2_parser>
	        {
	            typedef atom2_parser self_t;

	            const char* const firstChar;
	            mutable std::string label;

	            //! construction
	            atom2_parser(char * firstChar_) : firstChar(firstChar_), label("  ")
	            {}

	              //! tests whether or not the parsed character is a
	              //! valid and supported atom label
	              //! @param ch the parsed character to test
	              //! @return whether or not @a ch is a valid and supported
	              //! atom label
	            template <typename CharT>
	            bool test(CharT ch) const
	            {
	            	  // set label
					label[0] = *firstChar;
					label[1] = ch;
					  // check if two letter atom label is known
					return MoleculeUtil::getAtomData().find(label)
							!= MoleculeUtil::getAtomData().end();
	            }
	        };


	         /*!
	          * Dedicated parser for bond labels.
	          */
	        struct bond_parser : public NS_BOOSTSPIRIT::char_parser<bond_parser>
	        {
	            typedef bond_parser self_t;

	            //! construction
	            bond_parser() {}

	              //! tests whether or not the parsed character is a valid and
	              //! supported bond label
	              //! @param ch the parsed character to test
	              //! @return whether or not @a ch is a valid bond label
	            template <typename CharT>
	            bool test(CharT ch) const
	            {
	                return MoleculeUtil::getBondData().find(std::string(1,ch))
	                		!= MoleculeUtil::getBondData().end();
	            }
	        };

	        bond_parser const bond_p;



    	////////////// TYPEDEFs //////////////	
	    	
			
			  //! Utility helper class for parsing.
			class Atom_closure : public	NS_BOOSTSPIRIT::closure<	Atom_closure
								  							, std::string
								  							, std::string
								  							, int
								  							, int
					                                        , int >
			{
			public:
				typedef NS_BOOSTSPIRIT::closure<	Atom_closure
					, std::string
					, std::string
					, int
					, int
					, int > SuperClass;
					
				  //! bond label
				typename SuperClass::member1 blabel;
				  //! atom label
				typename SuperClass::member2 alabel;
				  //! atom index
				typename SuperClass::member3 cnt;
				  //! ring closure index
				typename SuperClass::member4 rc;
				  //! explicit H-atom count
				typename SuperClass::member5 hcnt;
			};

			  //! type of rule context of this class
			typedef NS_BOOSTSPIRIT::rule<	ScannerT
  											, typename Atom_closure::context_t
  									> rule_t;
  					
			// Typedefs for Local Data Structures
			typedef std::pair<int,int> bond_t;
			typedef std::vector<std::pair<bond_t,std::string> > bonds_t;
			struct AtomInfo {
				std::string label;
				int atomID;
				bool isAromatic;
				AtomInfo( const std::string label, const int atomID, const bool isAromatic )
				 : label(label), atomID(atomID), isAromatic(isAromatic)
				{}
			};
			typedef std::vector< AtomInfo > atoms_t;
			typedef std::vector<int> stack_t;
			typedef std::vector<std::pair<int,int> > hcount_t;
#if HAVE_UNORDERED_MAP > 0
			typedef std::unordered_map<int,int> rcs_t;
			typedef std::unordered_map<int,std::string> rcb_t;
#elif HAVE_TR1_UNORDERED_MAP > 0
			typedef std::tr1::unordered_map<int,int> rcs_t;
			typedef std::tr1::unordered_map<int,std::string> rcb_t;
#elif HAVE_GNU_HASH_MAP > 0
			class hash_int {
			public:

				size_t operator()(const int& v) const
				{
					return (size_t)v;
				}
				     
			};
			typedef __gnu_cxx::hash_map<int,int,hash_int> rcs_t;
			typedef __gnu_cxx::hash_map<int,std::string,hash_int> rcb_t;
#else
			typedef std::map<int,int> rcs_t;
			typedef std::map<int,std::string> rcb_t;
#endif


    	////////////// VARIABLES //////////////	
	    	
			  //! back reference to enclosing object for molecule creation
			SMILES_grammar const& self;
						
  			  // rules
			rule_t chain, simple_atom, complex_atom, ringclosure, bond, branch;
			
			NS_BOOSTSPIRIT::rule<ScannerT> smiles;
			NS_BOOSTSPIRIT::rule<ScannerT> simple_symbol, complex_symbol, bond_symbol, group_symbol;
			NS_BOOSTSPIRIT::rule<ScannerT> isotope, charge, chirality, hcount, name;
      
		     
		      // Local Data Structures
			int atom_count;
			  //! map of ring closure numbers to opening atom
			rcs_t rc_to_atom;
			  //! map of ring closure numbers to bond label of closing bond
			rcb_t rc_to_bond;
			atoms_t atoms;
			bonds_t bonds;
			stack_t stack;
			hcount_t explicitH;

	    	////////////// METHODS //////////////	

			  //! Function called by the parser to report the next atom.
			  //! @param atom the atom to add
			  //! @param alabel the label of the atom
			  //! @param blabel the label of the bond to the last reported atom
	          //! @throw std::invalid_argument in case a check fails
			void 
			memorize_atom( int atom, std::string& alabel, std::string& blabel )
				throw (std::invalid_argument);

	                  //! Function called by the parser to report the explicit hydrogen
	                  //! count for complex atoms
	                  //! @param atom the complex atom with explicit hydrogens
	                  //! @param hcount the number of explicit hydrogens
	                void
			memorize_explicit_H( int atom, int hcount);

			  //! Function called by the parser to report a ring closure.
			  //! @param rc the ring closure index that was closed
			  //! @param atom1 the first atom of the ring
	          //! @param blabel the optional given label of the ring closing bond
	          //! @throw std::invalid_argument in case a check fails
			void 
			memorize_rc( int rc, int atom1, std::string blabel )
				throw (std::invalid_argument);

			  //! Function called by the parser to report the opening of a 
			  //! molecule branching
			  //! @param atom the atom id where the branching occured
			void 
			open_branch( int atom );

			  //! Function called by the parser to report that the last
			  //! molecule branching is ended
			void 
			close_branch( void );
      
			  //! Resets the local data structures to allow for the parsing of 
			  //! the next SMILES string
			void 
			reset_data_structures( void );

	                  //! Function called by the parser that reports warnings
	                  //! to the standard error handle
	                void
			parser_warning( std::string msg );

	                  //! Function called by the parser to add explicit hydrogens
	                  //! from complex atoms to constructed graph
	                void
	                add_explicit_hydrogens(void);
      
     }; // struct definition
    

  }; // class SMILES_grammar


  
  } // namespace chem
} // namespace ggl

  // function implementations
#include "ggl/chem/SMILES_grammar.icc"


#endif /*SMILES_GRAMMAR_HH_*/
