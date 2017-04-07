#ifndef GGL_GRAPH_GML_GRAMMAR_HH_
#define GGL_GRAPH_GML_GRAMMAR_HH_

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
#error "GGL_GRAPH_GML_GRAMMAR : BOOST_SPIRIT_CLOSURE_LIMIT too low, has to be at least 5"
#endif

  // set phoenix limit if neccessary
#if !defined(PHOENIX_LIMIT)
#define PHOENIX_LIMIT 5
#elif PHOENIX_LIMIT < 5
#error "GGL_GRAPH_GML_GRAMMAR : PHOENIX_LIMIT too low, has to be at least 5"
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

#include "ggl/Graph.hh"

namespace ggl {


	/*! @brief Graph GML parser
	 *
	 * Parses a GML string representation of a ggl::Graph object.
	 * 
	 * Example :
	 * 
	 * \verbatim
	   ====== GRAPH TO ENCODE =======================

	     3(D) -1- 1(B)
	      |      /
	      2    3
	      |  /
	     2(C) -0- 0(A)

	   ======= GRAPH IN GML =========================

	   graph [
	     node [ id 0 label "A" ]
	     node [ id 1 label "B" ]
	     node [ id 2 label "C" ]
	     node [ id 3 label "D" ]
	     edge [ source 0 target 2 label "-0-" ]
	     edge [ source 1 target 3 label "-1-" ]
	     edge [ source 1 target 2 label "-3-" ]
	     edge [ source 2 target 3 label "-2-" ]
	   ]

	   ==============================================
	   \endverbatim
	 * 
	 * BNF grammar of GML (http://www.infosun.fim.uni-passau.de/Graphlet/GML/)
	 * 
	 * \verbatim
	    gml            ::= keyvalues
	    keyvalues      ::= keyvalue (keyvalue)
	    keyvalue       ::= key value
	    key            ::= ['a'-'z''A'-'Z']['a'-'z''A'-'Z''0'-'9']
	    value          ::= real | integer | string | list
	    real           ::= sign? digit  '.' digit+ mantissa
	    integer        ::= sign? digit+
	    string         ::= '"' instring '"'
	    list           ::= '[' keyvalues ']'
	    sign           ::= '+' | '-'
	    digit          ::= ['0'-'9']
	    mantissa       ::= ('E' | 'e') sign? digit+
	    instring       ::= ASCII-{'&', '"'} | '&' ['a'-'z''A'-'Z']  ';'
	   \endverbatim
	 *
	 * @author Christoph Flamm (c) 2008 http://www.tbi.univie.ac.at/~xtof/
	 * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */ 
	class Graph_GML_grammar
	  : public NS_BOOSTSPIRIT::grammar< Graph_GML_grammar >
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
		  LST_VAL  =  3  /* list value type    */
		};

		// typedefs
		typedef enum kv_values kv_values_t;

		typedef 
#if HAVE_UNORDERED_MAP > 0
		std::unordered_map<std::string, kv_values_t>
#elif HAVE_TR1_UNORDERED_MAP > 0
		std::tr1::unordered_map<std::string, kv_values_t>
#elif HAVE_GNU_HASH_MAP > 0
		__gnu_cxx::hash_map<std::string, kv_values_t, sgm::hash_string>
#else
		std::map<std::string, kv_values_t> 
#endif
		keys_map_t;

		typedef boost::tuple<int,int,std::string,std::string> edge_t;
		typedef boost::tuple<int,std::string,std::string> node_t;

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
		    _keys["id"]      = INT_VAL;
		    _keys["label"]   = STR_VAL;
		    _keys["source"]  = INT_VAL;
		    _keys["target"]  = INT_VAL;
		    _keys["graph"]   = LST_VAL;
		    _keys["node"]    = LST_VAL;
		    _keys["edge"]    = LST_VAL;
		    _keys["_node"]    = UNKNOWN;
		    _keys["_edge"]    = UNKNOWN;
		  }

		  kv_values_t lookup(std::string k) {
		    keys_map_t::iterator pos = _keys.find(k);

		    if (pos != _keys.end()) return (pos->second);
		    else return (UNKNOWN);
		  }

		  kv_values_t context(std::string k) {
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
	    Graph& g2fill;
	    
		  //! Access to the node label property_map of g2fill to set node labels
	    mutable 
	      boost::property_map< Graph, PropNodeLabel>::type
	    	nodeLabel;
	    
		  //! Access to the edge label property_map of g2fill to set edge labels
	    mutable 
	      boost::property_map< Graph, PropEdgeLabel>::type
	    	edgeLabel;
	    
		
	public:

		  //! Constructs the definitions of a GML graph grammar to parse
		  //! a GML graph string representation and to fill the encoded graph
		  //! into a given boost graph object.
		  //! @param toFill the graph object to add nodes and edges to
	    explicit Graph_GML_grammar( Graph & toFill );
	    
	      //! Parses a GML string and generates a Graph object
	      //! @param GML_string the string to parse
	      //! @return pair.first = the parsed graph object
	      //!         pair.second = -1 if parsing was successful,
	      //!         in error case it returns the string position that caused
	      //!         the parsing error
	    static
	    std::pair< Graph, int >
	    parseGraph( const std::string & GML_string );
	    

		  //! The definition of the GML grammar.
	    template <typename ScannerT>
	    struct definition
	    {
	    public:
	    	
	    	  //! Construction of the GML BNF grammar rules
	    	definition( Graph_GML_grammar const& self);
	    	
	    	  //! start parsing
	    	NS_BOOSTSPIRIT::rule<ScannerT> const&
	    	start() const;
	    	
	    protected:
	    	
	        // Variables explicitly initialized
	    	  //! back reference to enclosing object for molecule creation
	    	Graph_GML_grammar const& self;
	        typedef NS_BOOSTSPIRIT::rule<ScannerT
	        			, keyvalue_closure::context_t> keyvalue_t;
	        			
	        NS_BOOSTSPIRIT::rule<ScannerT> gml, keyvalues, value;
	        
	        keyvalue_t keyvalue, key, string, integer, real, list;
	        
	        // typedefs for local data structures

	        typedef std::vector<bool> boolstack_t;
	        typedef std::vector<std::string> keystack_t;
	        typedef std::set<edge_t,lt_edge> edges_t;
	        typedef std::set<node_t,lt_node> nodes_t;
	        
	        // helper functions

	        /***************************/
	        bool is_valid_node(node_t& a);
	        
	        /***************************/
	        bool is_valid_edge(edge_t& e);
	        
	        /***************************/
	        std::string spacer(int level);

	        /******************************/
	        void reset_data_structures(void);
	        
	        /****************************/
	        void clear_tmp_node_data(void);

	        /****************************/
	        void clear_tmp_edge_data(void);

	        /****************/
	        void dumpvec(void);

	        // semantic actions

	        /**************************/
	        void openList(std::string s);

	        /******************/
	        void closeList(void);

	        /*************/
	        void create_Graph(void);
	        
	        /*************/
	        void Dump(void);

	        /***************************************************************/
	        void dumpKeyValues(std::string k, std::string s, int i, double d);

	        /**********************/
	        void memorize_node(void);

	        /**********************/
	        void memorize_edge(void);

	        /***************************************************************/
	        void set_node_data(std::string k, std::string s, int i, double d);

	        /***************************************************************/
	        void set_edge_data(std::string k, std::string s, int i, double d);

	        /****************************************************************/
	        void keyValueAction(std::string k, std::string s, int i, double d);
	        
	        // local data structures

	        int level;
	        boolstack_t boolstack;
	        keystack_t keystack;
	        node_t tmp_node;
	        edge_t tmp_edge;
	        edges_t edges;
	        nodes_t nodes;
	        keys_map keys;
	     
	    	
	    };



	};


} // namespace ggl

// include implementation
#include "ggl/Graph_GML_grammar.icc"

#endif /*GRAPH_GML_GRAMMAR_HH_*/
