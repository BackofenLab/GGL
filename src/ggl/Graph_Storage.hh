#ifndef GGL_GRAPH_STORAGE_HH_
#define GGL_GRAPH_STORAGE_HH_

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>

#include "ggl/Graph.hh"

namespace ggl {



	  /*! @brief Interface graph storage
	   *
	   *  An abstract interface for a Graph object container for multiple
	   *  implementations. For instance, it is used by ggl::MR_ApplyRule to store the
	   *  graphs resulting from graph grammar ggl::Rule applications.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Graph_Storage {

	public:
		
		virtual ~Graph_Storage() 
		{}
		
		  //! Adds a given graph to the storage.
		  //! NOTE : The reported graphs might contain several independent
		  //! components! The Graph_Storage has to handle a separation if 
		  //! neccessary!
		  //! @param graph the Graph object to add.
		virtual
		void
		add( const Graph & graph ) = 0;
	};
	
} // namespace ggl


#endif /*GGL_GRAPH_STORAGE_HH_*/
