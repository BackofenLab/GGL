#ifndef GGL_GRAPH_GML_WRITER_HH_
#define GGL_GRAPH_GML_WRITER_HH_

#include <iostream>

#include "ggl/Graph.hh"


namespace ggl
{


  /*! @brief Graph GML writer
   *
   * Algorithm wrapper class to write a GGL graph in GGL GML format to stream.
   * 
   * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
   * 
   */
class Graph_GML_writer
{
public:
	  //! construction
	Graph_GML_writer();
	  //! destruction
	~Graph_GML_writer();
	
	
	  /*!
	   * Writes a boost graph in GGL GML format to stream. 
	   * 
	   * @param out the stream to write to
	   * @param graph the graph to print
	   * @param withSpaces if true newlines and whitespaces are used to write a
	   *        user friendly and readable output, otherwise a one-line output
	   *        with minimal space requirement is produced. 
	   * @param additionalGML if non-NULL, the according string is added to
	   *        the end of the graph GML encoding
	   */
	static 
	void
	write(	std::ostream& out
			, const Graph & graph
			, const bool withSpaces = true
			, const std::string* additionalGML = NULL );
	
};


} // namespace ggl


  // include implementation
#include "ggl/Graph_GML_writer.icc"

#endif /*GRAPH_GML_WRITER_HH_*/

