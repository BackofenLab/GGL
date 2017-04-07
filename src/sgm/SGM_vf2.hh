#ifndef SGM_VF2_H_
#define SGM_VF2_H_

#include "sgm/SubGraphMatching.hh"
#include "sgm/VF2_MatchingHandler.hh"

namespace sgm {


	  /*! @brief VF2-based subgraph matching
	   *
	   *  An sgm::SubGraphMatching interface implementation that uses the
	   *  implementation of the VF-2 algorithm for its computation.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class SGM_vf2 : public SubGraphMatching, public VF2_MatchingHandler {
		
	protected:
		
		typedef VF2_MatchingHandler::Label Label;
		typedef VF2_MatchingHandler::LabelComparator LabelComparator;
		typedef VF2_MatchingHandler::NodeData NodeData;
		typedef VF2_MatchingHandler::NodeDataVec NodeDataVec;
		typedef VF2_MatchingHandler::EdgeLabel EdgeLabel;
		typedef VF2_MatchingHandler::EdgeLabelVec EdgeLabelVec;
		typedef VF2_MatchingHandler::PatternMap PatternMap;
		typedef VF2_MatchingHandler::VF2_match_handler VF2_match_handler;
		typedef VF2_MatchingHandler::PatternVec PatternVec;
		typedef VF2_MatchingHandler::ReporterVec ReporterVec;

		  //! Comparator class for the node data used in the internal VF-2
		  //! data structures
		class NodeComparator: public LabelComparator {
		public:
			  //! the label defining the wildcard to use
			using LabelComparator::pWildcard;
			  //! access to the pattern currently matched
			const sgm::Pattern_Interface & pattern;
			  //! access to the target graph processed
			const sgm::Graph_Interface & target;
			  //! construction
			  //! @param pWildcard access to the wildcard label
			  //! @param pattern the pattern currently matched
			  //! @param target the target graph currently processed
			NodeComparator( const Label & pWildcard
						, const sgm::Pattern_Interface & pattern
						, const sgm::Graph_Interface & target);
			  //! comparison of the two pointers and objects.
			  //! @param pa first label pointer, assumed to be the pattern
			  //! @param pb second label pointer, assumed to be the target
			  //! @return true if the node data is compatible for matching
		    virtual bool compatible(void *pa, void *pb);
		};

		  //! Comparator class for the edge label used in the internal VF-2
		  //! data structures
		class EdgeLabelComparator: public LabelComparator {
		public:
			  //! the label defining the wildcard to use
			using LabelComparator::pWildcard;
			  //! construction
			  //! @param pWildcard access to the wildcard label
			EdgeLabelComparator( const Label & pWildcard );
			  //! comparison of the two pointers and objects.
			  //! @param pa first label pointer, assumed to be the pattern
			  //! @param pb second label pointer, assumed to be the target
			  //! @return true if pa and pb are Label pointers != NULL and pa
			  //! has a lexicographically smaller text 
		    virtual bool compatible(void *pa, void *pb);
		};


	public:


		////////////////////  SUB GRAPH MATCHING  //////////////////////////////


		  //! Performs exact sub graph matching to find maxHits occurrences of
		  //! the pattern graph within the target graph. Each hit is reported to
		  //! the Match_Reporter object.
		  //! @param pattern the pattern graph to search for
		  //! @param target the graph to search the pattern within
		  //! @param reporter all hits are reported to that object
		  //! @param maxHits the maximal number of hits to find
		  //! @return the number of exact matches found
		virtual
		size_t
		findMatches (	const Pattern_Interface & pattern,
						const Graph_Interface & target,
						Match_Reporter & reporter,
						const size_t maxHits );

		  //! Performs exact sub graph matching to find maxHits occurrences of
		  //! the pattern graphs within the target graph. Each hit is reported
		  //! to the Match_Reporter object.
		  //! NOTE, the first maxHits matches of the pattern graphs according
		  //! to their order in the patterns container are reported, i.e. first
		  //! all occurrences of the first pattern are identified. If this does
		  //! not exceed the maxHits limit, the next pattern is matched and so
		  //! on until either no pattern is left or the maxHits limit is
		  //! exceeded.
		  //! @param patterns the container of the pattern graphs to search for
		  //! @param target the graph to search the pattern within
		  //! @param reporter all hits are reported to that object
		  //! @param maxHits the maximal number of hits to find
		  //! @return the number of exact matches found
		virtual
		size_t
		findMatches (	const std::vector< const Pattern_Interface*> & patterns,
						const Graph_Interface & target,
						Match_Reporter & reporter,
						const size_t maxHits );

		  //! Performs exact sub graph matching to find maxHits occurrences of
		  //! the pattern graphs within the target graph. Each hit is reported
		  //! to the Match_Reporter object.
		  //! NOTE, the first maxHits matches of the pattern graphs according
		  //! to their order in the patterns container are reported, i.e. first
		  //! all occurrences of the first pattern are identified. If this does
		  //! not exceed the maxHits limit, the next pattern is matched and so
		  //! on until either no pattern is left or the maxHits limit is
		  //! exceeded.
		  //! @param patterns the container of the pattern graphs to search for
		  //! @param target the graph to search the pattern within
		  //! @param reporters each hit is reported to the corresponding object,
		  //!        the container has to have the same length as patterns
		  //! @param maxHits the maximal number of hits to find
		  //! @return the number of exact matches found
		virtual
		size_t
		findMatches (	const std::vector< const Pattern_Interface*> & patterns,
						const Graph_Interface & target,
						std::vector<Match_Reporter*> & reporters,
						const size_t maxHits );

	};

} // namespace sgm

#endif /*SGM_VF2_H_*/
