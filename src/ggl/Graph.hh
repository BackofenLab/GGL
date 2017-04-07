#ifndef GGL_GRAPH_HH_
#define GGL_GRAPH_HH_

#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>
#include <set>

#include "sgm/Graph_boost.hh"


  //! The ggl namespace contains classes needed to specify rules of a graph
  //! grammar and to apply these rules to graphs.
namespace ggl {


	  //! This boost graph property is used to determine the index of a given
	  //! node along the iterator order.
	typedef boost::vertex_index_t	PropNodeIndex;
	
	  //! Vector of node indices
	typedef std::vector< PropNodeIndex > NodeIndexVec;

	  //! Set of node indices
	typedef std::set< PropNodeIndex > NodeIndexSet;

	  //! This boost graph property is used to determine the label of a given
	  //! node.
	typedef boost::vertex_name_t	PropNodeLabel;

	  //! The properties available for the nodes of a Graph
	  //!
	  //! @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	typedef	boost::property<	PropNodeIndex, size_t
			, boost::property<	PropNodeLabel, std::string
				> >
					Graph_NodeProperties;
	
	  //! This boost graph property is used to determine the label of a given
	  //! edge.
	typedef boost::edge_name_t		PropEdgeLabel;

	  //! The properties available for the edges of a Graph
	  //!
	  //! @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	typedef	boost::property<	PropEdgeLabel, std::string 
				> 
					Graph_EdgeProperties;
	
	  //! The definition of undirected graphs that are handled by the graph
	  //! grammar rules of the GGL.
	  //!
	  //! @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  //!
	typedef boost::adjacency_list<
					boost::vecS,      				// store edges
					boost::vecS,       				// store vertices
					boost::undirectedS,				// is an undirected graph
					Graph_NodeProperties,  			// (atom symbols etc) 
					Graph_EdgeProperties   			// (edge symbols etc)
				> 
					Graph;


	  /*!
	   * Definition of a sgm::Graph_Interface wrapper around a ggl::Graph
	   * object.
	   *
	   * @author Martin Mann (c) 2012 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	typedef sgm::Graph_boost <	Graph \
							, PropNodeLabel \
							, PropEdgeLabel \
							, PropNodeIndex \
						> Graph_boost;

} // namespace ggl



 /*!
  * Writes a ggl::Graph object to stream in GML format
  * @param out the stream to write to
  * @param g the graph to be written in GML format
  * @return the changed input stream out
  */
std::ostream&
operator <<( std::ostream & out, const ggl::Graph& g );


#endif /*GRAPH_HH_*/
