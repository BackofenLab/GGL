
#ifndef GGL_CHEM_AP_DISABLED_HH_
#define GGL_CHEM_AP_DISABLED_HH_

#include <cassert>
#include <vector>
#include <set>
#include <algorithm>
#include <limits>


#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeUtil.hh"
#include "ggl/chem/AromaticityPerception.hh"

namespace ggl {
namespace chem {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


	 /*! @brief disables aromaticity prediction
	  *
	  * Dummy aromaticity perception implementation that actually performs NO
	  * AROMATICITY PERCEPTION at all and thus disables aromaticity perception
	  * wherever used.
	  *
	  * The current aromaticity assignment is preserved.
	  *
	  *  @author Martin Mann (c) 2014 http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class AP_disabled
		: public AromaticityPerception
	{
	public :

		  //! ring describing class
		typedef AromaticityPerception::RingDescriptor RingDescriptor;

		  //! Vector that holds ( index, weight ) pairs
		typedef std::vector< std::pair<size_t,double> > RingWeightVec;


	protected:

		  //! defines an edge; NOTE: use as ordered pairs, i.e. first <= second
		typedef AromaticityPerception::Edge Edge;
		  //! defines an aromatic ring via the set of aromatic edges
		typedef AromaticityPerception::EdgeSet EdgeSet;

		  //! the aromatic ring container to fill
		using AromaticityPerception::aromaticEdges;

		  //! container that will hold all rings of the current molecule
		  //! -> this list is later pruned to the rings of interest
		using AromaticityPerception::allRings;

	public:


		  /*!
		   * Construction
		   * @param model the aromaticity model to be used
		   */
		AP_disabled()
			: AromaticityPerception()
		{}

		  /*!
		   * Copy construction
		   * @param toCopy the object to make this a copy of
		   */
		AP_disabled( const AP_disabled & toCopy )
			: AromaticityPerception( toCopy )
		{}

		  /*!
		   * destruction
		   */
		virtual ~AP_disabled()
		{}


		  /*!
		   * Does nothing here...
		   *
		   * @param mol the molecule to check for rings
		   */
		virtual
		void
		findAllRings( const Molecule & mol )
		{
			// DO NOTHING
		}

		  /*!
		   * Push all currently aromatic edges to the
		   * aromaticEdges container.
		   *
		   * @param mol the molecule to check for aromatic rings
		   */
		void
		identifyAromaticEdges( const Molecule & mol )
		{
			boost::property_map<	Molecule , PropNodeIndex >
				::const_type nodeIndex = boost::get( PropNodeIndex(), mol );
			boost::property_map<	Molecule , PropEdgeLabel >
				::const_type edgeLabel = boost::get( PropEdgeLabel(), mol );

			  // get access to edge descriptors
			boost::graph_traits<Molecule>::edge_iterator ei, e_end;
			boost::tie(ei, e_end) = boost::edges( mol );

			  // check all edges if they have to be relabeled
			for (; ei != e_end; ++ei) {
				  // get bond data for this edge
				const MoleculeUtil::BondLabelData * bondData
							= MoleculeUtil::getBondData( edgeLabel[*ei] );
				assert( bondData != NULL /*unknown bond label*/);

				  // check if the bond is currently aromatic
				if ( bondData->isAromatic != 0 ) {
					  // generate edge descriptor
					Edge curEdge( nodeIndex[boost::source(*ei,mol)]
										, nodeIndex[boost::target(*ei,mol)] );
					if (curEdge.first > curEdge.second) {
						curEdge = Edge(curEdge.second, curEdge.first);
					}
					  // store to preserve bond aromaticity
					aromaticEdges.insert(curEdge);
				}
			}

		}

		 /*!
		  * Does nothing and thus fully disables aromaticity perception
		  *
		  * @param mol the molecule to correct (not changed)
		  * @param checkValence whether or not the valence of nodes and bonds
		  *        should be checked. (ignored)
		  *
		  * @throw std::runtime_error (never thrown)
		  */
		virtual
		void
		correctAromaticity( Molecule & mol, const bool checkValence )
			throw (std::runtime_error)
		{}


		 /*!
		  * Creates a heap copy of this instance that has to be deleted by the
		  * calling methods later on.
		  * @return a new instance that is a copy of this object
		  */
		virtual
		AP_disabled *
		clone() const
		{
			return new AP_disabled();
		}

	};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}} // namespace


#endif /* GGL_CHEM_AP_DISABLED_HH_ */
