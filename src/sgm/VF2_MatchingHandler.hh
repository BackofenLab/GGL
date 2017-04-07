
#ifndef SGM_VF2_MATCHINGHANDLER_HH_
#define SGM_VF2_MATCHINGHANDLER_HH_


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

#include <set>
#include <vector>

#include <vf2/argraph.h>
#include <vf2/argedit.h>

#include "sgm/Graph_Interface.hh"
#include "sgm/Pattern.hh"
#include "sgm/Match.hh"
#include "sgm/Match_Reporter.hh"

#include "sgm/MC_Node.hh"


namespace sgm {

	 /*! @brief General VF2-matching functionalities
	  *
	  * The general interface necessary to enable the use of the VF2 library
	  * for graph and sub graph matching. It covers the needed data structures
	  * and routines.
	  *
	  * @author Martin Mann (c) 2011 - http://www.bioinf.uni-freiburg.de/~mmann/
	  */
	class VF2_MatchingHandler {

	protected:

		  //! Data type that defines the the search mode of findMatches(..)
		  //! function
		enum ResultMode {FIND_ONE, FIND_ALL};


		  //! Label class used in the internal VF-2 data structures that just
		  //! indexes the known labels
		typedef int Label;

		  //! The data type that allows for the reuse of Label objects to be
		  //! more memory efficient.
		typedef
#if HAVE_UNORDERED_MAP > 0
			std::unordered_map< std::string, Label>
#elif HAVE_TR1_UNORDERED_MAP > 0
			std::tr1::unordered_map< std::string, Label>
#elif HAVE_GNU_HASH_MAP > 0
			__gnu_cxx::hash_map< std::string, Label, sgm::hash_string>
#else
			std::map< std::string, Label>
#endif
		 	 PatternMap;


		  //! Data handler for node information used for node comparison
		class NodeData {
		public:

			  //! container that holds node constraints to be verified
			typedef std::vector< const MC_Node* > NodeConstraints;

			  //! the node index of the node represented
			size_t id;
			  //! the label of the node
			Label label;
			  //! the out degree of the node
			size_t outDegree;
			  //! the number of self loops of the node
			size_t selfloops;
			  //! the node constraints to be checked
			NodeConstraints * nodeConstraints;
			  //! Default construction
			NodeData();
			  //! Construction
			  //! @param id the node index of the node represented
			  //! @param label the node label
			  //! @param outDegree the out degree of the node
			  //! @param selfloops the number of self loops of the node
			NodeData (	const size_t id
						, const Label label
						, const size_t outDegree
						, const size_t selfloops );
			  //! Copy construction
			  //! @param toCopy the object to make this a copy of
			NodeData( const VF2_MatchingHandler::NodeData & toCopy );
			  //! Destruction
			~NodeData();
		};

		  //! Comparator class for the node label used in the internal VF-2
		  //! data structures
		class LabelComparator: public vf2::AttrComparator {
		public:
			  //! the label defining the wildcard to use
			const Label pWildcard;
			  //! construction
			  //! @param pWildcard access to the wildcard label
			LabelComparator( const Label pWildcard );
			  //! comparison of the two pointers.
			  //! @param pa first label pointer
			  //! @param pb second label pointer
			  //! @return true if pa and pb are both != NULL and pa
			  //! equals pb
		    virtual bool compatible(void *pa, void *pb) = 0;
			  //! comparison of the two Label objects.
			  //! @param a first label
			  //! @param b second label
			  //! @return true if a and b are both Label != 0 and either
			  //! of them is a wildcard or both are equal
		    virtual bool compatibleLabel(const Label & a, const Label & b);
		};

		  //! Storage to enable fast memory clean up of NodeLabels.
		typedef std::vector< NodeData > NodeDataVec;

		  //! collector for parallel edge handling, since VF2 does not support
		  //! multiple edges between two nodes. Therefore, all labels of
		  //! parallel edges are compiled into one edge.
		typedef std::multiset<Label> EdgeLabel;

		  //! Storage to enable fast memory clean up of EdgeLabels.
		typedef std::vector< EdgeLabel* > EdgeLabelVec;


		//////////////////// VF2 MATCH HANDLING ///////////////////////////////

		  /*! @brief Handles match reporting from VF2-engine
		   *
		   * This class is needed to forward matches found by the VF-2 subgraph
		   * matching engine to a Match_Reporter instance.
		   *
		   */
		class VF2_match_handler {
		public:
			  //! The global Match_Reporter in use to write the found matchings to
			Match_Reporter& matchReporter;

			  //! The pattern reference used to write the found matchings
			  //! to matchReporter
			const Pattern_Interface& pattern;

			  //! The target reference used to write the found matchings
			  //! to matchReporter
			const Graph_Interface& target;

			  //! The maximal number of hits to find
			const size_t maxHitsToFind;

			  //! the number of matches for the current target to be matched
			size_t numOfMatches;

			  /*!
			   * Constructs a match handler for vf2::match(..) usage.
			   * @param matchReporter all matches will be reported to this
			   *        instance
			   * @param pattern the pattern graph needed to forward the report
			   *        to matchReporter
			   * @param target the target graph needed to forward the report
			   *        to matchReporter
			   * @param maxHitsToFind the maximal number of hits to find
			   */
			VF2_match_handler( Match_Reporter & matchReporter
								, const Pattern_Interface& pattern
								, const Graph_Interface& target
								, const size_t maxHitsToFind
								);

			  //! A static match reporting function used for the vf2-lib match(..)
			  //! call. It is reporting to the matchReporter member variable
			  //! access via x. The usage of the function will report
			  //! maxHitsToFind matches of the pattern in the target graph.
			  //! @param n length of the arrays q and t and so equal to the number
			  //!          of nodes of the pattern/query graph
			  //! @param q the matched indices from the pattern graph
			  //! @param t the node indices from the target graph where the indices
			  //!          from the pattern graph (q[i]) were matched to
			  //! @param x this will hold a pointer to a VF2_match_handler
			  //!          instance to access the reporter data
			  //! @return whether or not the search algorithm should search for
			  //!         another solution based on the number of matches seen
			  //!         (numOfMatches) and the targeted number (maxHitsToFind)
			static
			bool
			report_matches(int n, vf2::node_id *q, vf2::node_id *t, void *x);

		};


	public:

		VF2_MatchingHandler();
		virtual ~VF2_MatchingHandler();

	protected:

		  //! container of patterns to match
		typedef std::vector< const Pattern_Interface* > PatternVec;
		  //! container of match reporters to report matches to
		typedef std::vector< Match_Reporter* > ReporterVec;

		/*!
		 *  Performs exact sub-graph matching to find all occurences of the
		 *  pattern graph within the target graph. Each hit is reported to
		 *  the Match_Reporter object.
		 *  @param pattern the container of pattern graphs to search for
		 *  @param target the graph to search the pattern within
		 *  @param output container of match reporters, one for each pattern
		 *         to match, all hits of each pattern are reported to the
		 *         according output object
		 *  @param maxHits the maximal number of matches to find
		 *  @return the number of exact matches found
		 */
		template< class VF2STATE, class NODECOMPARE, class EDGECOMPARE >
		static
		size_t
		findMatches (	const PatternVec & pattern,
						const Graph_Interface & target,
						ReporterVec & output,
						const size_t maxHits );

		  /*!
		   *  Static function that is used by findMatches to fill the internally
		   *  used ARGEdit objects.
		   *  @param graph the graph to be 'filled' into the loader
		   *  @param loader the loader to fill
		   *  @param patternLabel the already known labels from the graph
		   *  @param nodeData the container that will hold the node data.
		   *  @param edgeLabels the container that will hold all edge labels.
		   *  @param isPattern If true, than graph will be used to fill
		   *         patternLabel and all nodes and edges in the loader will get
		   *         a label. If false, only known label from the pattern graph
		   *         are set in the loader, all other are set to NULL.
		   *  @return false if isPattern == false and at least one label in the
		   *         graph is not present in patternLabel, true otherwise
		   */
		static
		bool
		fillLoader(	const Graph_Interface& graph,
				  		vf2::ARGEdit & loader,
				  		PatternMap & patternLabel,
				  		NodeDataVec & nodeData,
				  		EdgeLabelVec & edgeLabels,
				  		bool isPattern );

		  /*!
		   * Checks whether or not a given pattern match fulfills the additional
		   * matching constraints requested by the pattern.
		   *
		   * @param pattern the pattern to be found
		   * @param target the graph the left side was found within, i.e. in
		   *   that the rule should be applied
		   * @param match contains the indices of the matched left side nodes in
		   *  the target graph. match[i] corresponds to the mapping of the i-th
		   *  vertex in the leftSide graph.
		   */
		static
		bool
		isApplicable (	const Pattern_Interface& pattern,
						const Graph_Interface & target,
						const Match & match );

		  /*!
		   * checks whether or not the pattern is compatible with the target
		   * based on the out degree information of the nodes.
		   *
		   * @param patternNodes the node data of the pattern
		   * @param targetNodes the node data of the target
		   * @return true, if the pattern node degree distribution is compatible
		   *        with the distribution of the target; false otherwise
		   */
		static
		bool
		areDegreeCompatible(	const NodeDataVec & patternNodes
								, const NodeDataVec & targetNodes );

	};

}

#include "sgm/VF2_MatchingHandler.icc"

#endif /* SGM_VF2_MATCHINGHANDLER_HH_ */
