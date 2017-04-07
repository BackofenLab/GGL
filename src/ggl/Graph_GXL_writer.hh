#ifndef GGL_GRAPH_GXL_WRITER_HH_
#define GGL_GRAPH_GXL_WRITER_HH_

#include <iostream>

#include "ggl/Graph.hh"


namespace ggl
{


  /*! @brief Graph GXL writer
   *
   * Algorithm wrapper class to write a boost graph in GXL format to stream.
   *
   *   http://www.gupro.de/GXL/
   * 
   * @author Martin Mann (c) 2011 http://www.bioinf.uni-freiburg.de/~mmann/
   * 
   */
class Graph_GXL_writer
{
public:
	  //! construction
	Graph_GXL_writer();
	  //! destruction
	~Graph_GXL_writer();
	
	
	  /*!
	   * Writes a boost graph in GXL format to stream.
	   * 
	   * @param out the stream to write to
	   * @param graph the graph to print
	   */
	static 
	void
	write(	std::ostream& out
			, const Graph & graph );
	
};


} // namespace ggl


  // include implementation
#include "ggl/Graph_GXL_writer.icc"

#endif /*GRAPH_GXL_WRITER_HH_*/

