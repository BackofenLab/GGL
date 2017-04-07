#ifndef GGL_CHEM_CLASSIDGRAPH_HH_
#define GGL_CHEM_CLASSIDGRAPH_HH_

#include "sgm/Graph_Interface.hh"

#include <sstream>

namespace ggl {
namespace chem {

	/**
	 * Graph wrapper that adds the atom node ID as class information to each
	 * atom node label. The class ID string length is the length of the node
	 * number string plus a given number of leading zeros. For instance a graph
	 * with 12 atoms nodes and two leading zeros will have class ID '0001' to
	 * '0012'.
	 *
	 * NOTE: existing class information is extended with the generated class ID
	 * string.
	 *
	 * @author Martin Mann - 2013
	 */
	class ClassIDGraph : public sgm::Graph_Interface {
	protected:

		const sgm::Graph_Interface & oriGraph;

		size_t classIdLength;

	public:

		 /**
		  * Construction
		  * @param molGraph the molecule graph to wrap
		  * @param leadingZeros the number of leading zeros in the class ID
		  */
		ClassIDGraph( const sgm::Graph_Interface & molGraph, const size_t leadingZeros )
		 :	oriGraph(molGraph)
			, classIdLength(0)
		{
			  // with positions for leading zeros
			classIdLength = 1 + leadingZeros;
			size_t divisor = 10;
			const size_t nodeNumber = oriGraph.getNodeNumber();
			while( nodeNumber / divisor > 0 ) {
				divisor *= 10;
				++classIdLength;
			}
		}

		/**
		 * Destruction
		 */
		virtual ~ClassIDGraph()
		{}


		  //! Access to the number of nodes of the graph
		  //! @return the overall node number
		virtual
		size_t
		getNodeNumber(void) const {
			return oriGraph.getNodeNumber();
		}

		  //! Access to iteration begin for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator to the first edge within the adjacency of i
		virtual
		OutEdge_iterator
		getOutEdgesBegin( const IndexType & i ) const {
			return oriGraph.getOutEdgesBegin(i);
		}

		  //! Access to iteration end for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator the end of the adjacency iteration of i
		virtual
		OutEdge_iterator
		getOutEdgesEnd( const IndexType & i ) const {
			return oriGraph.getOutEdgesEnd(i);
		}

		  //! Access to the label of a specified node
		  //! @param i the index of the node of interest
		  //! @return a string representation of the node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const {
			  // setup class ID generation
			std::ostringstream classID;
			classID.width(classIdLength);
			classID.fill('0');
			classID <<(i+1);

			  // get atom label
			const std::string oriLabel = oriGraph.getNodeLabel(i);
			  // check if class information already available
			if ( oriLabel.find(':') != std::string::npos ) {
				 // class ID present -> just append
				return oriLabel + classID.str();
			} else {
				 // no class ID present -> generate
				return oriLabel + std::string(":") + classID.str();
			}

		}

	};

}} // namespace

#endif /* GGL_CHEM_CLASSIDGRAPH_HH_ */
