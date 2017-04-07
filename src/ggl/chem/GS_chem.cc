
#include <vector>
#include <cassert>

#include <boost/graph/copy.hpp>

#include "ggl/chem/GS_chem.hh"

namespace ggl {
 namespace chem {

////////////////////////////////////////////////////////////////////////////////


	void
	GS_chem
	::add( const Molecule & graph )
	{
		if (boost::num_vertices(graph) == 0)
			return;

//std::cout <<"\n\n ##### GS_chem ####\n" <<std::endl;
		  // check graph for connected components
		std::vector<int> component(boost::num_vertices(graph));
		int num = boost::connected_components(graph, &component[0]);
//std::cout <<" GS_CHEM : num = "<<num <<std::endl;
//std::cout <<" GS_CHEM : ";
//for(size_t x=0; x<component.size();++x) {std::cout <<component.at(x) <<" ";}
//std::cout <<std::endl;
		assert( num > 0 /* no components found ... ?!? */);

		if (num == 1) {

//std::cout <<" GS_CHEM : just add" <<std::endl;
			  // forward reporting to sub class
			addMolecule( graph );
//std::cout <<" GS_CHEM : just add -> done" <<std::endl;

		} else {

//std::cout <<" GS_CHEM : split" <<std::endl;
			  // get index map
			const IndexMap idxMap = boost::get(PropNodeIndex(), graph);

//for (size_t i=0; i<component.size(); i++) {
//	std::cout <<" GS_CHEM : get(idxMap,"<<i<<") = "; std::cout.flush();
//	std::cout <<boost::get(idxMap,(Molecule::vertex_descriptor)i) <<std::endl;
//}

			  // create SMILES string for each component
			for (int comp=0; comp<num; ++comp) {

//std::cout <<" GS_CHEM : comp "<<comp <<std::endl;
				  // create subgraph of graph that represents current connected
				  // component

				  // set up filters for current connected component
				edge_is_in_component edgeFilter( comp, component, &idxMap, graph);
				node_is_in_component nodeFilter( comp, component, &idxMap);

				  // create boost graph representing the current component
				ComponentGraph compGraph(graph, edgeFilter, nodeFilter);

				  // the molecule graph representing the current component
				Molecule compMolGraph;

//std::cout <<" GS_CHEM : copy"<<std::endl;
				boost::copy_graph(compGraph, compMolGraph);
//std::cout <<" GS_CHEM : copy -> done"<<std::endl;

				  // forward reporting to sub class
				addMolecule( compMolGraph );

			}
//std::cout <<" GS_CHEM : split -> done" <<std::endl;
		}
	}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
 } // namespace chem
} //  namespace ggl
