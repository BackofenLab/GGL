#ifndef GGL_GS_STREAM_HH_
#define GGL_GS_STREAM_HH_

#include "ggl/Graph.hh"

#include "ggl/Graph_Storage.hh"

#include <iostream>

namespace ggl {



	  /*! @brief Graph storage within output stream
	   *
	   *  A ggl::Graph_Storage implementation that writes each added graph in string
	   *  representation to a given stream.
	   *
	   *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class GS_stream : public Graph_Storage
	{
		
	protected:
		
		std::ostream & out;
		
	public:
		
		  //! Construction
		  //! @param out the stream to write the graphs in text-form to
		GS_stream( std::ostream & out );
		
		virtual
		~GS_stream();
		
		
		  //! Writes a given graph to stream.
		  //! @param graph the Graph object to write.
		virtual
		void
		add( const Graph & graph );
	};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



} // namespace ggl

#include "ggl/GS_stream.icc"

#endif /*GS_STREAM_HH_*/
