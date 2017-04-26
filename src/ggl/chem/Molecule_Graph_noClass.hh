#ifndef GGL_CHEM_MOLECULE_GRAPH_NOCLASS_HH_
#define GGL_CHEM_MOLECULE_GRAPH_NOCLASS_HH_

#include "ggl/chem/Molecule.hh"

namespace ggl {
  namespace chem {

	  /*! @brief Wrapper for a molecule representing graph without atom class label
	   *
	   *  Represents a molecule class where class information is removed from
	   *  all atom labels
	   *
	   *  @author Martin Mann (c) 2017 http://www.bioinf.uni-freiburg.de/~mmann/
	   */
	class Molecule_Graph_noClass : public sgm::Graph_Interface {

	protected:
		
		//! wrapped molecules representation to strip the class label information from
		const ggl::chem::Molecule_Graph & molecules;
		
	public:

		//! Construction
		//! @param molecules the graph representing the molecules to wrap;
		//!    atom labels are reduced to their substrings up to the first
		//!    occurrence of the class label separator ':'
		Molecule_Graph_noClass( const ggl::chem::Molecule_Graph & molecules )
		 :	molecules(molecules)
		{}

		  //! Destruction
		virtual 
		~Molecule_Graph_noClass() {}
		
		  //! Access to the number of nodes of the graph
		  //! @return the overall node number 
		virtual
		size_t
		getNodeNumber(void) const
		{ return molecules.getNodeNumber(); }
		
		  //! Access to iteration begin for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator to the first edge within the adjacency of i
		virtual
		OutEdge_iterator
		getOutEdgesBegin( const IndexType & i ) const
		{ return molecules.getOutEdgesBegin( i ); }

		  //! Access to iteration end for the edge in the adjacency list of
		  //! a specified node
		  //! @param i the index of the node of interest
		  //! @return the iterator the end of the adjacency iteration of i
		virtual
		OutEdge_iterator
		getOutEdgesEnd( const IndexType & i ) const
		{ return molecules.getOutEdgesEnd( i ); }
		
		  //! Access to the label of a specified atom node where the original
		  //!    atom label is reduced to their substrings up to the first
		  //!    occurrence of the class label separator ':'
		  //! @param i the index of the node of interest
		  //! @return a string representation of the atom node label
		virtual
		std::string
		getNodeLabel(const IndexType & i) const
		{
			std::string origLabel = molecules.getNodeLabel(i);
			// return prefix of original atom label UP TO class label separator ':'
			return origLabel.substr( 0, origLabel.find(':') );
		}

		//! Provides access to the original graph without reduced atom labels.
		//! @return the original graph wrapped by this graph
		const ggl::chem::Molecule_Graph &
		getWithFullAtomLabels() const
		{ return molecules; }

	};

  } // namespace chem
} // namespace ggl



#endif /*GGL_CHEM_MOLECULE_GRAPH_NOCLASS_HH_*/
