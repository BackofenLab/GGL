
#include "ggl/Graph.hh"
#include "ggl/Graph_GML_writer.hh"

std::ostream&
operator <<( std::ostream & out, const ggl::Graph& g )
{

	  // convert graph into GML format
	ggl::Graph_GML_writer::write( out, g, true);

	return out;
}
