#ifndef GGL_RULE_GML_GRAMMAR_HH_
#define GGL_RULE_GML_GRAMMAR_HH_

#include <utility>
#include <vector>

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
#error "GGL_RULE_GML_GRAMMAR : BOOST_SPIRIT_CLOSURE_LIMIT too low, has to be at least 5"
#endif

  // set phoenix limit if neccessary
#if !defined(PHOENIX_LIMIT)
#define PHOENIX_LIMIT 5
#elif PHOENIX_LIMIT < 5
#error "GGL_RULE_GML_GRAMMAR : PHOENIX_LIMIT too low, has to be at least 5"
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

#include "sgm/MC_Node.hh"
#include "sgm/MC_Edge.hh"
#include "ggl/Rule.hh"
#include "ggl/Rule_GML_error.hh"

namespace ggl {



	/*! @brief Graph grammar rule parser
	 *
	 * Parses a GML string representation of a ggl::Rule object.
	 * 
	 * Example :
	 * 
	 * \verbatim
	   ====== RULE TO ENCODE =======================

	     3(D) -1- 1(B)               3(E) -1- 1(B)
	      |      /                    |        |
	      2    3           ==>        4        3
	      |  /                        |        |
	     2(C) -0- 0(A)               2(C)     4(D)

	   ======= RULE IN GML =========================

	   rule [
	           ruleID "Example rule"
	           context [
	                   node [ id 1 label "B" ]
	                   node [ id 2 label "C" ]
	                   edge [ source 1 target 3 label "-1-" ]
	           ]
	           left [
	                   node [ id 0 label "A" ]
	                   node [ id 3 label "D" ]
	                   edge [ source 0 target 2 label "-0-" ]
	                   edge [ source 1 target 2 label "-3-" ]
	                   edge [ source 2 target 3 label "-2-" ]
	           ]
	           right [
	                   node [ id 3 label "E" ]
	                   node [ id 4 label "D" ]
	                   edge [ source 2 target 3 label "-4-" ]
	                   edge [ source 1 target 4 label "-3-" ]
	           ]
	   ]

	   =============================================
	 * \endverbatim
	 * 
	 * BNF grammar of GML (http://www.infosun.fim.uni-passau.de/Graphlet/GML/)
	 * 
	   \verbatim
	    gml            ::= keyvalues
	    keyvalues      ::= keyvalue (keyvalue)
	    keyvalue       ::= key value
	    key            ::= ['a'-'z''A'-'Z']['a'-'z''A'-'Z''0'-'9']
	    value          ::= real | integer | string | operator | list
	    real           ::= sign? digit  '.' digit+ mantissa
	    integer        ::= sign? digit+
	    string         ::= '"' instring '"'
	    operator       ::= '>' | '<' | '=' | '!'
	    list           ::= '[' keyvalues ']'
	    sign           ::= '+' | '-'
	    digit          ::= ['0'-'9']
	    mantissa       ::= ('E' | 'e') sign? digit+
	    instring       ::= ASCII-{'&', '"'} | '&' ['a'-'z''A'-'Z']  ';'
	 * \endverbatim
	 *
	 * @author Christoph Flamm (c) 2008 http://www.tbi.univie.ac.at/~xtof/
	 * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */ 
	class Rule_GML_grammar
	  : public NS_BOOSTSPIRIT::grammar< Rule_GML_grammar >
	{
	protected:
		
		  //! enummeration of value types in gml key-value pairs
		enum kv_values {
		  Y_CONTXT = -3, /* works as outer key context */
		  N_CONTXT = -2, /* ignore as outer context key */
		  UNKNOWN  = -1, /* unknown value type */
		  INT_VAL  =  0, /* integer value type */
		  STR_VAL  =  1, /* string value type  */
		  DBL_VAL  =  2, /* double value type  */
		  LST_VAL  =  3,  /* list value type    */
		  STRLST_VAL  =  4  /* string list value type    */
		};

		// typedefs
		typedef enum kv_values kv_values_t;
		typedef boost::tuple<int,int,std::string,std::string> edge_t;
		typedef boost::tuple<int,std::string,std::string> node_t;
#if HAVE_UNORDERED_MAP > 0
		typedef std::unordered_map<std::string, kv_values_t> keys_map_t;
#elif HAVE_TR1_UNORDERED_MAP > 0
		typedef std::tr1::unordered_map<std::string, kv_values_t> keys_map_t;
#elif HAVE_GNU_HASH_MAP > 0
		typedef __gnu_cxx::hash_map<std::string, kv_values_t, sgm::hash_string> keys_map_t;
#else
		typedef std::map<std::string, kv_values_t> keys_map_t;
#endif

		struct lt_edge : public std::binary_function<edge_t, edge_t, bool> {
		  bool operator()(edge_t a, edge_t b) {
		    if (a.get<0>() != b.get<0>()) return ( a.get<0>() < b.get<0>() );
		    if (a.get<1>() != b.get<1>()) return ( a.get<1>() < b.get<1>() );
		    if (a.get<2>() != b.get<2>()) return ( a.get<2>() < b.get<2>() );
		    if (a.get<3>() != b.get<3>()) return ( a.get<3>() < b.get<3>() );

		    return (false);
		  }
		};

		struct lt_node : public std::binary_function<node_t, node_t, bool> {
		  bool operator()(node_t a, node_t b) {
		    if (a.get<0>() != b.get<0>()) return ( a.get<0>() < b.get<0>() );
		    if (a.get<1>() != b.get<1>()) return ( a.get<1>() < b.get<1>() );
		    if (a.get<2>() != b.get<2>()) return ( a.get<2>() < b.get<2>() );

		    return (false);
		  }
		};
		
		struct keyvalue_closure : NS_BOOSTSPIRIT::closure<keyvalue_closure,
		    std::string, std::string, int, double>
		{
			typedef NS_BOOSTSPIRIT::closure<	keyvalue_closure
								, std::string
								, std::string
								, int
								, double > SuperClass;
								
			SuperClass::member1 key;
			SuperClass::member2 strval;
			SuperClass::member3 numval;
			SuperClass::member4 dblval;
		};
		
		class keys_map {
		  public:

		  keys_map() { inizialize(); }

		  void inizialize(void) {
		    _keys["label"]    = STR_VAL;
		    _keys["ruleID"]   = STR_VAL;
		    _keys["wildcard"] = STR_VAL;
		    _keys["id"]       = INT_VAL;
		    _keys["source"]   = INT_VAL;
		    _keys["target"]   = INT_VAL;
		    _keys["graph"]    = LST_VAL;
		    _keys["rule"]     = LST_VAL;
		    _keys["context"]  = LST_VAL;
		    _keys["left"]     = LST_VAL;
		    _keys["right"]    = LST_VAL;
		    _keys["node"]     = LST_VAL;
		    _keys["edge"]     = LST_VAL;
		    _keys["_node"]    = UNKNOWN;
		    _keys["_edge"]    = UNKNOWN;
		    _keys["constrainAdj"]	= LST_VAL;
		    _keys["constrainNode"]	= LST_VAL;
		    _keys["constrainNoEdge"]	= LST_VAL;
		    _keys["constrainEdge"]	= LST_VAL;
		    _keys["count"]			= INT_VAL;
		    _keys["op"]				= STR_VAL;
		    _keys["nodeLabels"]		= LST_VAL;
		    _keys["edgeLabels"]		= LST_VAL;
		    _keys["copyAndPaste"]	= LST_VAL;
		  }

		  kv_values_t lookup( const std::string & k) {
		    keys_map_t::iterator pos = _keys.find(k);

		    if (pos != _keys.end()) return (pos->second);
		    else return (UNKNOWN);
		  }

		  kv_values_t context( const std::string & k) {
		    keys_map_t::iterator pos = _keys.find("_"+k);

		    if (pos != _keys.end()) return (N_CONTXT);
		    else return (Y_CONTXT);
		  }

		  private:
		  keys_map_t _keys;
		};


	protected:
		
		  //! The boost core graph object that is filled to represent the next 
		  //! parsed Rule.
	    Rule::CoreGraph& g2fill;
	    
	      //! The rule id that is filled.
	    std::string& ruleID;

	      //! The additional constraints to be fulfilled by the rule to be filled
	    std::vector< sgm::Pattern_Interface::Match_Constraint* > & ruleConstraints;

	      //! The copy-and-Paste operation container to be filled
	    Rule::CopyAndPasteOperations & copyAndPaste;

	      //! The wildcard label that is filled if present.
	    std::string& wildcard;

		  //! Access to the node label property_map of g2fill to set node labels
	    mutable 
	      boost::property_map<Rule::CoreGraph
	      					, Rule::NodeLabelProperty>::type
	    	nodeLabel;
	    
		  //! Access to the right node label property_map of g2fill to set
	      //! changing node labels
	    mutable 
	      boost::property_map<Rule::CoreGraph
	      					, Rule::NodeRightLabelProperty>::type
	    	nodeLabelRight;
	    
		  //! Access to the node Rule context property_map of g2fill 
	    mutable 
	      boost::property_map<Rule::CoreGraph
	      					, Rule::NodeContextProperty>::type
	    	nodeContext;
	    
		  //! Access to the edge label property_map of g2fill to set edge labels
	    mutable 
	      boost::property_map<Rule::CoreGraph
	      					, Rule::EdgeLabelProperty>::type
	    	edgeLabel;
	    
		  //! Access to the edge Rule context property_map of g2fill 
	    mutable 
	      boost::property_map<Rule::CoreGraph
	      					, Rule::EdgeContextProperty>::type
	    	edgeContext;
	    
		
	public:

		  //! Constructs the definitions of a GML graph grammar to parse
		  //! a GML graph string representation and to fill the encoded graph
		  //! into a given boost graph object.
		  //! @param toFill the Rule core graph object to add nodes and edges to
		  //! @param ruleID the Rule ID to set
	      //! @param ruleConstraints the additional constraints that
	      //!        have to be met by the rule
	      //! @param copyAndPaste the copy-and-Paste operation container to fill
	      //! @param wildcard the wildcard string to be filled. NOTE: if no
	      //!        wildcard was parsed this string is not changed, so check
	      //!        for a change to know if a wildcard was parsed or not!
	    explicit Rule_GML_grammar(	Rule::CoreGraph& toFill
	    							, std::string& ruleID
	    							, std::vector< sgm::Pattern_Interface::Match_Constraint* > & ruleConstraints
	    							, Rule::CopyAndPasteOperations & copyAndPaste
	    							, std::string& wildcard );
	    
	      //! Parses a GML string and generates a Rule object
	      //! @param GML_string the string to parse
	      //! @return pair.first = the graph encoding of the molecule
	      //!         pair.second = -1 if parsing was successfull,  
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
	      //! @throw ggl::Rule_GML_error if parsing errors occur
	    static
	    std::pair< Rule, int >
	    parseRule( const std::string & GML_string ) throw (Rule_GML_error);
	    

		  //! The definition of the GML grammar.
	    template <typename ScannerT>
	    struct definition
	    {
	    public:
	    	
	    	  //! Construction of the GML BNF grammar rules
	    	definition( Rule_GML_grammar const& self);
	    	
	    	  //! start parsing
	    	NS_BOOSTSPIRIT::rule<ScannerT> const&
	    	start() const;
	    	
	    protected:
	    	
	        // Variables explicitly initialized
	    	  //! back reference to enclosing object for molecule creation
	    	Rule_GML_grammar const& self;
	        typedef NS_BOOSTSPIRIT::rule<ScannerT
	        			, keyvalue_closure::context_t> keyvalue_t;
	        			
	        NS_BOOSTSPIRIT::rule<ScannerT> gml, keyvalues, value;
	        
	        keyvalue_t keyvalue, key, string, integer, real, relop, strlist, list;
	        
	        // typedefs for local data structures

	        typedef std::vector<bool> boolstack_t;
	        typedef std::vector<std::string> keystack_t;
	        typedef std::multiset<edge_t,lt_edge> edges_t;
	        typedef std::set<node_t,lt_node> nodes_t;
	        
	        // helper functions

	        /***************************/
	        bool is_valid_node(const node_t& a);

	        /***************************/
	        bool is_valid_edge(const edge_t& e);
	        
	        /***************************/
	        bool is_valid_MC_NodeAdjacency(const sgm::MC_NodeAdjacency& c);
	        
	        /***************************/
	        bool is_valid_MC_NodeLabel(const sgm::MC_NodeLabel& c);

	        /***************************/
	        bool is_valid_MC_NoEdge(const sgm::MC_NoEdge& c);

	        /***************************/
	        bool is_valid_MC_EdgeLabel(const sgm::MC_EdgeLabel& c);

	        /***************************/
	        bool is_valid_copyAndPaste(const Rule::RuleCnP& cnp);

	        /***************************/
	        std::string spacer(int level);

	        /******************************/
	        void reset_data_structures(void);
	        
	        /****************************/
	        void clear_tmp_node_data(void);

	        /****************************/
	        void clear_tmp_edge_data(void);

	        /****************************/
	        void clear_tmp_MC_NodeAdjacency(void);

	        /****************************/
	        void clear_tmp_MC_NodeLabel(void);

	        /****************************/
	        void clear_tmp_MC_NoEdge(void);

	        /****************************/
	        void clear_tmp_MC_EdgeLabel(void);

	        /****************************/
	        void clear_tmp_copyAndPaste(void);

	        /****************/
	        void dumpvec(void);

	        // semantic actions

	        /**************************/
	        void openList(std::string s);

	        /******************/
	        void closeList(void);

	        /*************/
	        void create_Rule(void);
	        
	        /*************/
	        void Dump(void);

	        /***************************************************************/
	        void dumpKeyValues(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_node_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_edge_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_NodeAdjacency_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_NodeAdjacencyNL_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_NodeAdjacencyEL_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_NodeLabel_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_NodeLabelNL_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_NoEdge_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_EdgeLabel_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_MC_EdgeLabelEL_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_copyAndPaste_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_copyAndPaste_EL_data(std::string k, std::string s, int i, double d);

	        /****************************************************************/
	        void keyValueAction(std::string k, std::string s, int i, double d);
	        
	        // local data structures

	        int level;
	        boolstack_t boolstack;
	        keystack_t keystack;
	        node_t tmp_node;
	        edge_t tmp_edge;
	        std::string tmp_ruleID;
	        std::string tmp_wildcard;
	        edges_t edges;
	        nodes_t nodes;
	        keys_map keys;

		      //! The current sgm::MC_NodeAdjacency rule constraint to be filled
		    sgm::MC_NodeAdjacency tmp_MC_NodeAdjacency;
	     
		      //! The current sgm::MC_NodeLabel rule constraint to be filled
		    sgm::MC_NodeLabel tmp_MC_NodeLabel;

		      //! The current sgm::MC_NoEdge rule constraint to be filled
		    sgm::MC_NoEdge tmp_MC_NoEdge;

		      //! The current sgm::MC_EdgeLabel rule constraint to be filled
		    sgm::MC_EdgeLabel tmp_MC_EdgeLabel;

		      //! The current copy-and-paste operation to be filled
		    Rule::RuleCnP tmp_copyAndPaste;

	    	
	    };

	protected:
		


		//! Trims leading and tailing quotes from a string
		//! @param s the string to trim
		//! @return the trimmed string
	  static
	  std::string
	  trimQuotes( const std::string& s );


	};


} // namespace ggl

// include implementation
#include "ggl/Rule_GML_grammar.icc"

#endif /*RULE_GML_GRAMMAR_HH_*/
