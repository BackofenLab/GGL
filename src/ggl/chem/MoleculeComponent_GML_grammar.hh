#ifndef GGL_CHEM_MOLECULECOMPONENT_GML_GRAMMAR_HH_
#define GGL_CHEM_MOLECULECOMPONENT_GML_GRAMMAR_HH_

#include <utility>
#include <vector>
#include <string>
#include <stdexcept>

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

  // set spirit closure limit if neccessary
#if !defined(BOOST_SPIRIT_CLOSURE_LIMIT)
#define BOOST_SPIRIT_CLOSURE_LIMIT 5
#elif BOOST_SPIRIT_CLOSURE_LIMIT < 5
#error "GGL_CHEM_MOLECULECOMPONENT_GML_GRAMMAR : BOOST_SPIRIT_CLOSURE_LIMIT too low, has to be at least 5"
#endif

  // set phoenix limit if neccessary
#if !defined(PHOENIX_LIMIT)
#define PHOENIX_LIMIT 5
#elif PHOENIX_LIMIT < 5
#error "GGL_CHEM_MOLECULECOMPONENT_GML_GRAMMAR : PHOENIX_LIMIT too low, has to be at least 5"
#endif

#include <boost/version.hpp>
#if BOOST_VERSION >= 103800
#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/phoenix1.hpp>
#include <boost/spirit/include/classic_actor.hpp>
#define NS_BOOSTSPIRIT boost::spirit::classic
#else
#include <boost/spirit.hpp>
#include <boost/spirit/phoenix.hpp>
#include <boost/spirit/actor.hpp>
#define NS_BOOSTSPIRIT boost::spirit
#endif

#include "ggl/chem/MoleculeComponent.hh"
#include "sgm/Pattern.hh"
#include "ggl/chem/MC_MC_Node.hh"

namespace ggl {
namespace chem {


	/*! @brief MoleculeComponent parser
	 *
	 * Parses a GML string representation of a
	 * ggl::MoleculeDecomposition::MoleculeComponent object. This includes its
	 * properties as well as the additional constraints needed for matching.
	 * 
	 * Example :
	 * 
	 * \verbatim
	   ======= GRAPH IN GML =========================

	   molcomp [
	     description " '-Cl' (attached to a primary carbon with no other clorine atoms attached)"
	     priority 4
	     energy  -11.7
	     node [ id 0 label "C" ]
	     node [ id 1 label "Cl" ]
	     edge [ source 0 target 1 label "-" ]
	     constrainAdj [
		    id 0
		    op =
		    count 1
		    nodeLabels [ label "Cl" ]
	     ]
	     constrainAdj [
	        id 0
	        op =
	        count 2
	        nodeLabels [ label "H" ]
	     ]
	   ]

	   ==============================================
	   \endverbatim
	 * 
	 * @author Martin Mann (c) 2010 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */ 
	class MoleculeComponent_GML_grammar
	  : public NS_BOOSTSPIRIT::grammar< MoleculeComponent_GML_grammar >
	{
	protected:
		

		  //! type for mapping integers to size_t
		typedef 
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map<int, size_t>
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map<int, size_t>
#else
		std::map<int, size_t>
#endif
		MapIntSizeT;
		

	protected:
		
		  //! The boost core graph object that is filled to represent the next 
		  //! parsed Rule.
	    MoleculeComponent& toFill;

	public:
	    
		  //! Constructs the definitions of a GML graph grammar to parse
		  //! a GML graph string representation and to fill the encoded graph
		  //! into a given boost graph object.
		  //! @param toFill the object to add nodes and edges to
	    explicit MoleculeComponent_GML_grammar( MoleculeComponent & toFill );
	    
	      //! Parses a GML string and generates a MoleculeComponent::PatternGraph object
	      //! @param GML_string the string to parse
	      //! @return pair.first = the graph encoding of the molecule
	      //!         pair.second = -1 if parsing was successfull,
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
	      //! @throw std::invalid_argument in case a check fails
	    static
	    std::pair< MoleculeComponent, int >
	    parseGML( const std::string & GML_string ) throw (std::invalid_argument);
	    

		  //! The definition of the GML grammar.
	    template <typename ScannerT>
	    struct definition
	    {
	    public:
	    	
	    	  //! Construction of the GML BNF grammar rules to parse a
	    	  //! MoleculeComponent
	    	  //! @param self access to the calling grammar
	    	definition( const MoleculeComponent_GML_grammar & self );
	    	
	    	  //! start parsing
	    	NS_BOOSTSPIRIT::rule<ScannerT> const&
	    	start() const;
	    	
	    protected:
	    	

	          //! the molecule component to be filled
	        MoleculeComponent& toFill;

			  //! Access to the node label property_map of toFill to set node labels
			boost::property_map< MoleculeComponent::PatternGraph, PropNodeLabel>::type
		    	nodeLabel;

			  //! Access to the edge label property_map of g2fill to set edge labels
			boost::property_map< MoleculeComponent::PatternGraph, PropEdgeLabel>::type
		    	edgeLabel;

	          //! the mapping of node IDs in the GML notation and their
	          //! corresponding node IDs in the created pattern graph
	        MapIntSizeT nodeMapping;
	        			

	          //! the rules to be parsed
	        NS_BOOSTSPIRIT::rule<ScannerT> molcomp, content, node, edge,
					compIDs, constrainAdjacency, constrainLabel, ringFragment;


	        // temporary data structures
			int curNodeID, curEdgeFromID, curEdgeToID;
			std::string curNodeLabel, curEdgeLabel;
			std::string curRingFragmentTypeString;
			// temporary constraint objects
			MC_MC_NodeAdjacency constrAdj;
			char constrAdjOp, constrLabelOp;
			MC_MC_NodeLabel constrLabel;
			MoleculeComponent::RingFragmentList curRingFragmentList;

			enum WhatList {
				Fill_constrAdj_NL,
				Fill_constrAdj_EL,
				Fill_constrLabel_NL,
				Fill_compIDs,
				Fill_ringFragments
			};

	        // helper functions

	        /*!
	          * Resets toFill to enable the detection of missing information
	          * within final_checks()
	          */
	        void clear_toFill(void);

	        /*!
	          * Performs final checks on toFill to ensure that the parse was
	          * correct
	          *
	          * @throw std::invalid_argument in case a check fails
	          */
	        void final_checks(void) throw (std::invalid_argument);

	        /*!
	          * Stores a node of the MoleculeComponent pattern graph.
	          * @throw std::invalid_argument in case a check fails
	          */
	        void store_node(void) throw (std::invalid_argument);

	        /*!
	          * Stores an edge of the MoleculeComponent pattern graph.
	          * @throw std::invalid_argument in case a check fails
	          */
	        void store_edge(void) throw (std::invalid_argument);

	        /*!
	          * Stores a node adjacency constraint of the MoleculeComponent
	          */
	        void store_constrAdj(void);

	        /*!
	          * Stores a node label constraint of the MoleculeComponent
	          */
	        void store_constrLabel(void);

	        /*!
	          * Inserts a string into a given string set
	          * @param list the encoding which list to be fill
	          */
	        void insert_to_list( WhatList list );

	    }; // end of description


	}; // end of MoleculeComponent_GML_grammar


} // namespace chem
} // namespace ggl

// include implementation
#include "ggl/chem/MoleculeComponent_GML_grammar.icc"

#endif /*GGL_CHEM_MOLECULECOMPONENT_GML_GRAMMAR_HH_*/
