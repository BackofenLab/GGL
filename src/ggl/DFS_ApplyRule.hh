
#ifndef GGL_DFS_APPLYRULE_HH_
#define GGL_DFS_APPLYRULE_HH_

#include "sgm/SGM_vf2.hh"

#include "ggl/Graph.hh"
#include "ggl/Graph_Storage.hh"
#include "ggl/Rule.hh"

namespace ggl {

///////////////////////////////////////////////////////////////////////////////



	 /*! Depths-First-Search for graph grammars
	  *
	  * A generic Depths-First-Search (DFS) implementation that performs a
	  * recursive graph grammar rule application, i.e. it starts a new
	  * iteration of rule matching and application
	  * for each reported result graph of the last iteration.
	  *
	  * The recursion end is to be defined by an instance of the RecursionEnd.
	  *
	  *  @author Martin Mann (c) 2011 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class DFS_ApplyRule : public Graph_Storage {

	public:

		  /*!
		   * Handler to decide if a solution was found, a trace back is to be
		   * performed or another iteration of DFS is needed.
		   */
		class DFS_Visitor {

		public:

			  //! The possible decisions returned by the visitor.
			enum Decision { SOLUTION_STOP, SOLUTION_TRACEBACK, FAILURE_STOP, FAILURE_TRACEBACK, CONTINUE };

			  //! Destruction
			virtual
			~DFS_Visitor() {}

			  /*!
			   * Decides on the status of the current DFS, i.e. if the search
			   * is to be continued, a solution was found, or a traceback has
			   * to be done.
			   *
			   * @param graph the current graph to be decided on
			   * @return the status of the DFS
			   */
			virtual
			Decision
			status( const Graph & graph ) = 0;

		};

		 //! Container that stores a DFS trace.
		typedef std::vector< const Graph* > SearchTrace;

	public:

		 /*! Construction
		  */
		DFS_ApplyRule( );

		 //! Destruction
		virtual ~DFS_ApplyRule();

		 /*!
		  * Performs a DFS when applying the given rule onto the defined start
		  * graph.
		  *
		  * @param rules the rules to apply
		  * @param startGraph the start graph for the DFS
		  * @param solutionStorage the container where to report the solution to
		  * @param visitor the visitor instance that guides the DFS. Each result
		  *        graph is checked with the visitor to decide if a solution was
		  *        found, a trace back is needed or further DFS is to be done.
		  * @param searchTrace if non-NULL, DFS will copy the trace of the DFS
		  *        search to the container, i.e. the sequence of graphs
		  *        generated along the search
		  * @param doSymmBreak whether or not symmetry breaking should be done
		  * @param sgm the sub graph matcher to be used; if NULL an sgm::SGM_vf2
		  *        object is used per default
		  *
		  * @return whether or not a solution was found
		  */
		bool
		findSolution( const std::vector< Rule > & rules
					, const Graph & startGraph
					, Graph_Storage & solutionStorage
					, DFS_Visitor & visitor
					, SearchTrace * searchTrace = NULL
					, const bool doSymmBreak = true
					, sgm::SubGraphMatching * sgm = NULL );


	public:  // GRAPH_STORAGE INTERFACE -> starts next DFS iteration


		  /*!
		   * The Graph_Storage interface is used to trigger a new DFS iteration.
		   * Thus, each added graph is checked if the DFS is to be aborted or
		   * not. In case the search is to be extended, another rule application
		   * iteration is started.
		   *
		   * @param graph the graph that might start another DFS iteration or
		   * is either a solution or DFS dead end.
		   */
		void
		add( const Graph & graph );


	protected:

		  //! the rule pattern currently applied
		std::vector< const sgm::Pattern_Interface * > rulePattern;
		  //! the solution storage currently used
		Graph_Storage * solutionStorage;
		  //! the DFS visitor currently used
		DFS_Visitor * visitor;
		  //! the sub graph matcher currently used
		sgm::SubGraphMatching * sgm;
		  //! the match reporter to be informed about each match
		std::vector< sgm::Match_Reporter * > reporter;
		  //! the final trace to be filled if non-NULL
		SearchTrace * finalTrace;
		  //! the current trace maintained by the current search
		SearchTrace currentTrace;
		  //! whether or not a solution was found during this search
		bool solutionFound;

		  //! Dummy class to trigger DFS abortion via exception handling.
		class DFS_ABORTION {};


	};

///////////////////////////////////////////////////////////////////////////////

}


#include "ggl/DFS_ApplyRule.icc"


#endif /* GGL_DFS_APPLYRULE_HH_ */
