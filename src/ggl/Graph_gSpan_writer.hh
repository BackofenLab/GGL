#ifndef GGL_GRAPH_GSPAN_WRITER_HH_
#define GGL_GRAPH_GSPAN_WRITER_HH_

#include <iostream>

#include "ggl/Graph.hh"


namespace ggl
{


  /*! @brief Graph gSpan writer
   *
   * Algorithm wrapper class to write a boost graph in gSpan format to stream.
   * 
   * e.g.
   *
   *  ####################
   *  t
   *  v 1 A
   *  v 2 B
   *  e 1 2 lala
   *  ####################
   *
   * @author Martin Mann (c) 2013 http://www.bioinf.uni-freiburg.de/~mmann/
   * 
   */
class Graph_gSpan_writer
{
public:
	  //! construction
	Graph_gSpan_writer()
	{}

	  //! destruction
	~Graph_gSpan_writer()
	{}
	
	
	  /*!
	   * Writes a boost graph in GGL GML format to stream. 
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



#endif /*GGL_GRAPH_GSPAN_WRITER_HH_*/

