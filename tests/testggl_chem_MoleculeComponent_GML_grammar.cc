
#include <iostream>
#include <sstream>

#include "ggl/chem/MoleculeComponent.hh"
#include "ggl/chem/MoleculeComponent_GML_grammar.hh"

#include "sgm/Graph_boost.hh"

#include "utilPrintGraph_Interface.icc"

using namespace ggl;
using namespace chem;

const std::string gml[] = {
"  molcomp [\n"
"     description \" '-Cl' (attached to a primary carbon with no other chlorine atoms attached)\"\n"
"     priority 4\n"
"     energy  -11.7\n"
"     node [ id 0 label \"*\" ]\n"
"     node [ id 1 label \"Cl\" ]\n"
"     edge [ source 0 target 1 label \"-\" ]\n"
"     compIDs [ id 1 ]\n"
"     constrainAdj [\n"
"       id 0\n"
"       op =\n"
"       count 1\n"
"       nodeLabels [ label \"Cl\" ]\n"
"     ]\n"
"     constrainAdj [\n"
"       id 0\n"
"       op =\n"
"       count 2\n"
"       nodeLabels [ label \"H\" ]\n"
"     ]\n"
"     constrainNode [\n"
"       id 0\n"
"       op =\n"
"       nodeLabels [ label \"C\" label \"c\" ]\n"
"     ]\n"
"  ]"
	,
"  molcomp [\n"
"     description \" '-S-' (participating in a ring)\"\n"
"     priority 2\n"
"     energy  0.72\n"
"     node [ id 0 label \"*\" ]\n"
"     node [ id 1 label \"S\" ]\n"
"     node [ id 2 label \"*\" ]\n"
"     edge [ source 0 target 1 label \"-\" ]\n"
"     edge [ source 1 target 2 label \"-\" ]\n"
"     compIDs [ id 1 ]\n"
"     ringFragment [ id 0 id 1 id 2 type nonaromatic ]\n"
"  ]"
	,
	"" // marks the end of the array
};



void print( const MoleculeComponent & m ) {
	std::cout
			<<" + description = '" <<m.description <<"'" <<std::endl
			<<" + priority    = '" <<m.priority <<"'" <<std::endl
			<<" + energy      = '" <<m.freeEnergy <<"'" <<std::endl
			<<" + pattern     =\n\n" <<m.pattern <<std::endl
			<<" + compIDs     = ";
	for (MoleculeComponent::NodeSet::const_iterator n = m.compIDs.begin(); n != m.compIDs.end(); n++) {
		std::cout <<*n <<" ";
	}
	std::cout
			<<std::endl
			<<" + constraints = " <<m.constraints.size() <<std::endl
			;
	for (size_t i=0; i<m.constraints.size(); i++) {
		if (dynamic_cast<sgm::MC_NodeAdjacency*>(m.constraints[i]) != NULL) {
			sgm::MC_NodeAdjacency* tmp = dynamic_cast<sgm::MC_NodeAdjacency*>(m.constraints[i]);
			std::cout <<"  - cAdj "
						<<tmp->constrainedNodeID <<" ";
			switch(tmp->op) {
			case sgm::MC_NodeAdjacency::MC_EQ : std::cout <<"="; break;
			case sgm::MC_NodeAdjacency::MC_L : std::cout <<"<"; break;
			case sgm::MC_NodeAdjacency::MC_G : std::cout <<">"; break;
			default: std::cout <<tmp->op; break;
			}
			std::cout <<" " <<tmp->count <<" NL={ ";
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator it= tmp->nodeLabels.begin(); it != tmp->nodeLabels.end(); it++) {
				std::cout <<(*it) <<" ";
			}
			std::cout <<"} EL={";
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator it= tmp->edgeLabels.begin(); it != tmp->edgeLabels.end(); it++) {
				std::cout <<(*it) <<" ";
			}
			std::cout <<"}\n";
		}
		if (dynamic_cast<sgm::MC_NodeLabel*>(m.constraints[i]) != NULL) {
			sgm::MC_NodeLabel* tmp = dynamic_cast<sgm::MC_NodeLabel*>(m.constraints[i]);
			std::cout <<"  - cNod "
						<<tmp->constrainedNodeID <<" ";
			std::cout <<" NL={ ";
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator it= tmp->nodeLabels.begin(); it != tmp->nodeLabels.end(); it++) {
				std::cout <<(*it) <<" ";
			}
			std::cout <<"}\n";
		}
	}
	std::cout
			<<" + rings       = " <<m.ringFragments.size() <<std::endl;
	for (MoleculeComponent::RingFragVec::const_iterator r = m.ringFragments.begin(); r != m.ringFragments.end(); r++) {
		std::cout <<"  - ";
		for (MoleculeComponent::RingFragmentList::const_iterator n = r->fragment.begin(); n!=r->fragment.end(); n++) {
			std::cout <<*n <<" ";
		}
		switch( r->type ) {
		case MoleculeComponent::RF_aromaticHydrocarbon : std::cout <<" = aromatic hydrocarbon"; break;
		case MoleculeComponent::RF_heteroaromatic : std::cout <<" = heteroaromatic"; break;
		case MoleculeComponent::RF_nonaromatic : std::cout <<" = non-aromatic"; break;
		case MoleculeComponent::RF_undefined : std::cout <<" = type undefined"; break;
		}
		std::cout <<std::endl;
	}
}

int main() {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"    ggl::chem::MoleculeComponent_GML_grammar  \n"
				<<"==============================================\n" 
				<<std::endl;


	for (size_t i = 0; !(gml[i].empty()); i++) {
		std::cout	<<"\n######### next GML #########\n"
					<<gml[i]
					<<"\n############################\n"
					<<std::endl;

		std::cout <<"\n-->MoleculeComponent_GML_grammar::parse( GML )" <<std::endl;

		std::pair<MoleculeComponent, int> ret = MoleculeComponent_GML_grammar::parseGML(gml[i]);

		if (ret.second < 0 )
		{

			std::cout <<"\n result MoleculeComponent : " <<std::endl;
			print(ret.first);
			std::cout <<"\n" <<std::endl;
		} else {
			std::cout <<"\n PARSING ERROR : at input position " <<ret.second
					<<" = '" <<gml[i].at(ret.second) <<"' "
					<<" within \n'"
					<<gml[i].substr(std::max(0,((int)ret.second)-10),20)
					<<"'\n"
					<<std::endl;
		}
	}



	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<std::endl;
	return 0;
}
