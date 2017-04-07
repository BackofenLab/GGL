#ifndef GGL_GS_STL_HH_
#define GGL_GS_STL_HH_

#include "ggl/Graph.hh"

#include "ggl/Graph_Storage.hh"

#include <vector>

namespace ggl {

#define GGL_GS_STL_PUSHALLT_TEMPLATE \
	template <	class STL_CONTAINER >

#define GGL_GS_STL_PUSHALLT_TYPE \
	GS_STL_pushAllT <	STL_CONTAINER >

	/*! @brief Graph storage within STL push-back container
	 *
	 *  A ggl::Graph_Storage implementation that pushs each added graph
	 *  to the end of a provided STL container using its "push_back" method.
	 *
	 *  @tparam STL_CONTAINER a push-back STL container for Graph objects.
	 *
	 *  @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	template <	class STL_CONTAINER = typename std::vector< Graph > >
	class GS_STL_pushAllT : public Graph_Storage
	{
	public:
		
		  //! type of STL container of graphs this graph storage uses
		typedef STL_CONTAINER	storage_type;
		
	protected:
		
		  //! the STL container to write the graphs to
		STL_CONTAINER & storage;
		
	public:
		
		  //! Construction
		  //! @param storage the STL container to write the graphs to
		GS_STL_pushAllT( STL_CONTAINER & storage );
		
		  //! destruction
		virtual
		~GS_STL_pushAllT();
		
		
		  //! Adds the given graph to the internal STL container.
		  //! @param graph the Graph object to add.
		virtual
		void
		add( const Graph & graph );
	};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	  //! Wrapper typedef around the default parameters of GS_STL_pushAllT<..>
	typedef GS_STL_pushAllT<> GS_STL_pushAll;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ggl


namespace ggl {

#define GGL_GS_STL_PUSHALLPT_TEMPLATE \
	template <	class STL_CONTAINER >

#define GGL_GS_STL_PUSHALLPT_TYPE \
	GS_STL_pushAllPT <	STL_CONTAINER >


	/*!
	 * @brief Stores new graph pointers in STL push-back container.
	 *
	 * A ggl::Graph_Storage implementation that pushs a 'new' allocated
	 * copy of each added graph
	 * to the end of a provided STL container using its "push_back" method.
	 *
	 * @tparam STL_CONTAINER a push-back STL container for Graph pointer.
	 *
	 * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	template <	class STL_CONTAINER = typename std::vector< Graph* > >
	class GS_STL_pushAllPT : public Graph_Storage
	{
	public:
		
		  //! type of STL container of graphs this graph storage uses
		typedef STL_CONTAINER	storage_type;
		
	protected:
		
		STL_CONTAINER & storage;
		
	public:
		
		  //! Construction
		  //! @param storage the STL container to write the graphs to
		GS_STL_pushAllPT( STL_CONTAINER & storage );
		
		  //! destruction
		virtual
		~GS_STL_pushAllPT();
		
		
		  //! Adds the given graph to the internal STL container.
		  //! @param graph the Graph object to add.
		virtual
		void
		add( const Graph & graph );
	};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	  //! Wrapper typedef around the default parameters of GS_STL_pushAllPT<..>
	typedef GS_STL_pushAllPT<> GS_STL_pushAllP;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ggl


#include "sgm/GraphMatching.hh"

namespace ggl {

#define GGL_GS_STL_PUSHUNIQUET_TEMPLATE \
	template <	class STL_CONTAINER >

#define GGL_GS_STL_PUSHUNIQUET_TYPE \
	GS_STL_pushUniqueT <	STL_CONTAINER >


	/*!
	 * @brief Unique storage of graphs in STL push-back container
	 *
	 * A ggl::Graph_Storage implementation that pushs the added graph
	 * to the end of a provided STL container !!! ONLY if it is not only
	 * present in the container !!!, using its "push_back" method.
	 *
	 * The check is done using the provided sgm::GraphMatching object.
	 * Thus the more elements are in the container the more tests have to be
	 * done which results in slower runtimes!
	 *
	 * @tparam STL_CONTAINER a push-back STL container for Graph objects.
	 *
	 * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	template <	class STL_CONTAINER = typename std::vector< Graph > >
	class GS_STL_pushUniqueT : public Graph_Storage
	{
	public:
		
		  //! type of STL container of graphs this graph storage uses
		typedef STL_CONTAINER	storage_type;
		
	protected:
		
		  //! the STL container to write the graphs to
		STL_CONTAINER & storage;
		
		  //! the graph matcher used to decide if the added
		  //! graph is already present in the storage
		sgm::GraphMatching& matcher;
		
		  //! the wildcard to use for the matching
		const std::string * wildcard;

	public:
		
		  //! Construction
		  //! @param storage the STL container to write the graphs to
		  //! @param matcher the graph matcher used to decide if the added 
		  //!        graph is already present in the storage
		GS_STL_pushUniqueT(	STL_CONTAINER & storage
							, sgm::GraphMatching& matcher );

		  //! Construction
		  //! @param storage the STL container to write the graphs to
		  //! @param matcher the graph matcher used to decide if the added
		  //!        graph is already present in the storage
		  //! @param wildcard the wildcard to use for the matching.
		GS_STL_pushUniqueT(	STL_CONTAINER & storage
							, sgm::GraphMatching& matcher
							, const std::string& wildcard );

		  //! Copy construction
		  //! @param toCopy the object to make this a copy of
		GS_STL_pushUniqueT(	const GGL_GS_STL_PUSHUNIQUET_TYPE &toCopy );
		
		  //! destruction
		virtual
		~GS_STL_pushUniqueT();
		
		
		  //! Adds the given graph to the internal STL container.
		  //! @param graph the Graph object to add.
		virtual
		void
		add( const Graph & graph );
	};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	  //! Wrapper typedef around the default parameters of GS_STL_pushUniqueT<..>
	typedef GS_STL_pushUniqueT<> GS_STL_pushUnique;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ggl


#include "sgm/GraphMatching.hh"

namespace ggl {

#define GGL_GS_STL_PUSHUNIQUEPT_TEMPLATE \
	template <	class STL_CONTAINER >

#define GGL_GS_STL_PUSHUNIQUEPT_TYPE \
	GS_STL_pushUniquePT <	STL_CONTAINER >


	/*!
	 * @brief Unique storage of new graph pointers in STL push-back container
	 *
	 * A ggl::Graph_Storage implementation that pushs the a 'new' allocated
	 * copy of the added graph
	 * to the end of a provided STL container !!! ONLY if it is not only
	 * present in the container !!!, using its "push_back" method.
	 *
	 * The check is done using the provided sgm::GraphMatching object.
	 * Thus the more elements are in the container the more tests have to be
	 * done which results in slower runtimes!
	 *
	 * @tparam STL_CONTAINER a push-back STL container for Graph pointer.
	 *
	 * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	 */
	template <	class STL_CONTAINER = typename std::vector< Graph* > >
	class GS_STL_pushUniquePT : public Graph_Storage
	{
	public:
		
		  //! type of STL container of graphs this graph storage uses
		typedef STL_CONTAINER	storage_type;
		
	protected:
		
		  //! the STL container to write the graphs to
		STL_CONTAINER & storage;
		
		  //! the graph matcher used to decide if the added
		  //! graph is already present in the storage
		sgm::GraphMatching& matcher;
		
		  //! the wildcard to use for the matching
		const std::string * wildcard;

	public:
		
		  //! Construction
		  //! @param storage the STL container to write the graphs to
		  //! @param matcher the graph matcher used to decide if the added 
		  //!        graph is already present in the storage
		GS_STL_pushUniquePT(	STL_CONTAINER & storage
							, sgm::GraphMatching& matcher);
		
		  //! Construction
		  //! @param storage the STL container to write the graphs to
		  //! @param matcher the graph matcher used to decide if the added
		  //!        graph is already present in the storage
		  //! @param wildcard the wildcard to use for the matching.
		GS_STL_pushUniquePT(	STL_CONTAINER & storage
							, sgm::GraphMatching& matcher
							, const std::string& wildcard );

		  //! Copy construction
		  //! @param toCopy the object to make this a copy of
		GS_STL_pushUniquePT(	const GGL_GS_STL_PUSHUNIQUEPT_TYPE &toCopy );

		  //! destruction
		virtual
		~GS_STL_pushUniquePT();
		
		
		  //! Adds the given graph to the internal STL container.
		  //! @param graph the Graph object to add.
		virtual
		void
		add( const Graph & graph );
	};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	  //! Wrapper typedef around the default parameters of GS_STL_pushUniquePT<..>
	typedef GS_STL_pushUniquePT<> GS_STL_pushUniqueP;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ggl

  // include implementation
#include "ggl/GS_STL.icc"

#endif /*GS_STL_HH_*/
