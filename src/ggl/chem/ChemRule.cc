
#include "ggl/chem/ChemRule.hh"

namespace ggl {
  namespace chem {
  
  
////////////////////////////////////////////////////////////////////////////////
//	CONSTANTS
////////////////////////////////////////////////////////////////////////////////


	const size_t 
	ChemRule_Constants::C_AtomLabelInvalid = 5;

	const size_t 
	ChemRule_Constants::C_BondLabelInvalid = 7;

	const size_t 
	ChemRule_Constants::C_UnbalancedElectrons = 11;

	const size_t
	ChemRule_Constants::C_AtomComplexWithH = 13;

	const size_t
	ChemRule_Constants::C_BondLoop = 17;

	const size_t
	ChemRule_Constants::C_EdgeMultiplicity = 19;

	const size_t
	ChemRule_Constants::C_NodeInsertDeletion = 23;


////////////////////////////////////////////////////////////////////////////////



	ChemRule &
	ChemRule
	::operator =( const ChemRule & toCopy )
	{
		  // call assignment routine of super class
		static_cast<SuperClass*>(this)->operator=( toCopy );
		  // remove old transition state if present
		if (transitionState != NULL)
			delete transitionState;
		  // copy transition state if present, otherwise set to NULL
		transitionState = (toCopy.transitionState == NULL)
							? NULL
							: new TransitionState(*(toCopy.transitionState));
		  // return *this access
		return *this;
	}



////////////////////////////////////////////////////////////////////////////////



	size_t
	ChemRule
	:: isConsistent( void ) const
	{
		  // run consistency checks of superclass
		size_t retCode = SuperClass::isConsistent();

		  // check if no node is created or destroyed
		if (!checkNodeInDel( core )) {
			retCode *= C_NodeInsertDeletion;
		}

		  // check if all bonds are correctly defined
		if (!checkEdgeInDel( core )) {
			retCode *= C_EdgeMultiplicity;
		}

		  // check if all node label are ok
		if (!checkAtomLabel( core )) {
			retCode *= C_AtomLabelInvalid;
		} else {

			////  CHECKS THAT REQUIRE CORRECT ATOM LABEL

			  // check if electron change is balanced
			if (!checkElectronChange( core )) {
				retCode *= C_UnbalancedElectrons;
			}

			  // check if an atom label is complex and contains addition proton
			  // information -> currently not allowed, all hydrogen atom have to
			  // be explicit
			if (!checkAtomComplexWithH( core )) {
				retCode *= C_AtomComplexWithH;
			}

		}

		  // check if all bond label are ok
		if (!checkBondLabel( core )) {
			retCode *= C_BondLabelInvalid;
		}

		  // check if bond ring present
		if (!checkBondLoop( core )) {
			retCode *= C_BondLoop;
		}

		  // return final consistency status
		return retCode;
	}


////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	:: decodeConsistencyStatus(	const size_t consistencyCode
								, std::ostream& errorStream
								, const bool completeCheck )
	{
		if (consistencyCode == C_Consistent) {
			return true;
		}
		  // temporary error code handling
		size_t errorCode = consistencyCode;
		errorStream	<<"The chemical rule is not consistent because :\n";

		  // get error code handling of super class
		SuperClass::decodeConsistencyStatus( consistencyCode, errorStream, false );

		if (errorCode % C_NoID == 0) {
			errorStream <<" + contains no ruleID information \n";
			errorCode /= C_NoID;
		}
		if (errorCode % C_NodeInsertDeletion == 0) {
			errorStream <<" + contains node insertion/deletion \n";
			errorCode /= C_NodeInsertDeletion;
		}
		if (errorCode % C_EdgeMultiplicity == 0) {
			errorStream <<" + contains too many edges for one of the bonds \n";
			errorCode /= C_EdgeMultiplicity;
		}
		if (errorCode % C_AtomLabelInvalid == 0) {
			errorStream <<" + contains SMILES incompatible atom labels \n";
			errorCode /= C_AtomLabelInvalid;
		}
		if (errorCode % C_BondLabelInvalid == 0) {
			errorStream <<" + contains SMILES incompatible bond labels \n";
			errorCode /= C_BondLabelInvalid;
		}
		if (errorCode % C_UnbalancedElectrons == 0) {
			errorStream <<" + contains unbalanced electron changes \n";
			errorCode /= C_UnbalancedElectrons;
		}
		if (errorCode % C_AtomComplexWithH == 0) {
			errorStream <<" + contains implicit protons (H) within complex atom label, currently not supported \n";
			errorCode /= C_AtomComplexWithH;
		}
		if (errorCode % C_BondLoop == 0) {
			errorStream <<" + an atom forms a bond with itself \n";
			errorCode /= C_BondLoop;
		}
		if (completeCheck && errorCode != 1) {
			errorStream <<" + undescribed error code "<<errorCode<<" \n";
		}

		return false;
	}


////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkNodeInDel( const ChemRule::CoreGraph & coreGraph )
	{
		  // check that all nodes are context --> no 'atom creation'
		CoreGraph::vertex_iterator vIt, vItEnd;
		boost::property_map<	CoreGraph,
										NodeContextProperty >::const_type
			nodeLevel = boost::get( NodeContextProperty(), coreGraph );
		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(coreGraph); vIt != vItEnd; ++vIt) {
			  // check if current node is in context or label is changed
			if (	nodeLevel[*vIt] != RULE_CONTEXT
					&& nodeLevel[*vIt] != RULE_LABEL_CHANGE )
			{
				  // current node not part of context, i.e. created or destroyed
				return false;
			}
		}
		  // all nodes are in context
		return true;
	}



////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkEdgeInDel( const ChemRule::CoreGraph & coreGraph )
	{
		const size_t LEFT = 1, RIGHT = 2, COMPLETE = 3;
		std::map< std::pair< size_t, size_t >, size_t > edges;

		  // property access
		boost::property_map<	CoreGraph,
										EdgeContextProperty >::const_type
			edgeContext = boost::get( EdgeContextProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeIndexProperty >::const_type
			nodeIndex = boost::get( NodeIndexProperty(), coreGraph );
		  // check all bonds
		CoreGraph::edge_iterator eIt, eItEnd;
		for (boost::tie(eIt,eItEnd) = boost::edges(coreGraph); eIt != eItEnd; ++eIt) {
			size_t from = nodeIndex[ boost::source(*eIt,coreGraph) ];
			size_t to = nodeIndex[ boost::target(*eIt,coreGraph) ];
			std::pair< size_t, size_t > curEdge(from<to?from:to,from<to?to:from);

			  // check if context
			switch( edgeContext[*eIt] ) {
			case RULE_CONTEXT :
			{
				  // if context no other edges allowed
				if( edges.find(curEdge) != edges.end() )
					return false;
			} break;
			  // check if left
			case RULE_LEFT_SIDE :
			{
				  // check if known
				if (edges.find(curEdge) == edges.end()) {
					edges[curEdge] = LEFT;
				} else {
					  // check if second edge was a right edge --> label change
					if (edges[curEdge] == RIGHT) {
						  // no further edges allowed
						edges[curEdge] = COMPLETE;
					} else {
						  // a non-right edge was also present
						return false;
					}
				}
			} break;
			  // check if right
			case RULE_RIGHT_SIDE :
			{
				  // check if known
				if (edges.find(curEdge) == edges.end()) {
					edges[curEdge] = RIGHT;
				} else {
					  // check if second edge was a left edge --> label change
					if (edges[curEdge] == LEFT) {
						  // no further edges allowed
						edges[curEdge] = COMPLETE;
					} else {
						  // a non-right edge was also present
						return false;
					}
				}
			} break;
			default:
				assert(false /* edge with RULE_LABEL_CHANGE context --> not allowed */);
				break;
			}
		}
		  // all edge indels are correct
		return true;
	}



////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkAtomLabel( const ChemRule::CoreGraph & coreGraph )
	{
		CoreGraph::vertex_iterator vIt, vItEnd;
		boost::property_map<	CoreGraph,
										NodeContextProperty >::const_type
			nodeLevel = boost::get( NodeContextProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeLabelProperty >::const_type
			leftLabel = boost::get( NodeLabelProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeRightLabelProperty >::const_type
			rightLabel = boost::get( NodeRightLabelProperty(), coreGraph );

		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(coreGraph); vIt != vItEnd; ++vIt) {
			  // check left side / context label
			if (!MoleculeUtil::isValidAtomLabel( leftLabel[*vIt] )) {
				return false;
			}
			  // check if label change --> go to next handling
			if ( nodeLevel[*vIt] == RULE_LABEL_CHANGE
					&& !MoleculeUtil::isValidAtomLabel( rightLabel[*vIt] ))
			{
				return false;
			}
		}
		  // all nodes are ok
		return true;
	}

////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkAtomComplexWithH( const ChemRule::CoreGraph & coreGraph )
	{
		CoreGraph::vertex_iterator vIt, vItEnd;
		boost::property_map<	CoreGraph,
										NodeContextProperty >::const_type
			nodeLevel = boost::get( NodeContextProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeLabelProperty >::const_type
			leftLabel = boost::get( NodeLabelProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeRightLabelProperty >::const_type
			rightLabel = boost::get( NodeRightLabelProperty(), coreGraph );

		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(coreGraph); vIt != vItEnd; ++vIt) {
			  // access shortcut
			const std::string* curNodeLabel = &(leftLabel[*vIt]);
			  // check if the left/contetx label shows implicit protons
			if (	MoleculeUtil::getAtom(*curNodeLabel).compare("H") != 0
					&& MoleculeUtil::getProtons(*curNodeLabel) > 0 )
			{
				return false;
			}
			  // check if label change --> go to next handling
			if ( nodeLevel[*vIt] == RULE_LABEL_CHANGE ) {
				curNodeLabel = &(rightLabel[*vIt]);
				  // check if the right label shows implicit protons
				if (	MoleculeUtil::getAtom(*curNodeLabel).compare("H") != 0
						&& MoleculeUtil::getProtons(*curNodeLabel) > 0 )
				{
					return false;
				}
			}
		}
		  // all nodes are ok
		return true;
	}

////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkBondLabel( const ChemRule::CoreGraph & coreGraph )
	{
		CoreGraph::edge_iterator eIt, eItEnd;
		boost::property_map<	CoreGraph,
										EdgeLabelProperty >::const_type
			edgeLabel = boost::get( EdgeLabelProperty(), coreGraph );

		  // check all bonds
		for (boost::tie(eIt,eItEnd) = boost::edges(coreGraph); eIt != eItEnd; ++eIt) {
			  // check if current node is in context or label is changed
			  // check left side / context label
			if (!MoleculeUtil::isValidBondLabel( edgeLabel[*eIt] )) {
				return false;
			}
		}
		  // all bonds are ok
		return true;
	}


////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkBondLoop( const ChemRule::CoreGraph & coreGraph )
	{
		CoreGraph::edge_iterator eIt, eItEnd;

		  // check all bonds
		for (boost::tie(eIt,eItEnd) = boost::edges(coreGraph); eIt != eItEnd; ++eIt) {
			  // check if source and target are identical, i.e. ring bond found
			if ( boost::source(*eIt, coreGraph) == boost::target(*eIt, coreGraph)) {
				return false;
			}
		}
		  // all bonds are ok
		return true;
	}


////////////////////////////////////////////////////////////////////////////////



	bool
	ChemRule
	::checkElectronChange( const ChemRule::CoreGraph & coreGraph )
	{
		CoreGraph::vertex_iterator vIt, vItEnd;
		CoreGraph::out_edge_iterator aIt, aItEnd;
		CoreGraph::edge_iterator eIt, eItEnd;

		boost::property_map<	CoreGraph,
										NodeContextProperty >::const_type
			nLevel = boost::get( NodeContextProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeLabelProperty >::const_type
			nLeftLabel = boost::get( NodeLabelProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										NodeRightLabelProperty >::const_type
			nRightLabel = boost::get( NodeRightLabelProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										EdgeContextProperty >::const_type
			eLevel = boost::get( EdgeContextProperty(), coreGraph );
		boost::property_map<	CoreGraph,
										EdgeLabelProperty >::const_type
			edgeLabel = boost::get( EdgeLabelProperty(), coreGraph );

		  // check all atom nodes
		for (	boost::tie(vIt,vItEnd)
				= boost::vertices(coreGraph); vIt != vItEnd; ++vIt)
		{
			int inc = 0;
			int dec = 0;
			int chargeChange = 0;
			int protonChange = 0;

			bool isAromaticNode = false;

			  // get label change based information
			switch (nLevel[*vIt]) {
			case RULE_LABEL_CHANGE : {
				chargeChange = MoleculeUtil::getCharge(nRightLabel[*vIt])
								- MoleculeUtil::getCharge(nLeftLabel[*vIt]);
				  // handle special case of charge changes for protons
				if (MoleculeUtil::getAtom(nLeftLabel[*vIt])=="H") {
					chargeChange *= -1;
				}
				protonChange = MoleculeUtil::getProtons(nRightLabel[*vIt])
								- MoleculeUtil::getProtons(nLeftLabel[*vIt]);
				  // check whether left or right atom label is an aromatic one
				isAromaticNode = MoleculeUtil::getAtomData( nLeftLabel[*vIt] )->isAromatic == 1
								|| MoleculeUtil::getAtomData( nRightLabel[*vIt] )->isAromatic == 1;
				break;
			}
			case RULE_CONTEXT :
			case RULE_LEFT_SIDE :
			case RULE_RIGHT_SIDE :
			{
				  // check whether atom label is an aromatic one
				isAromaticNode = MoleculeUtil::getAtomData( nLeftLabel[*vIt] )->isAromatic == 1;
				break;
			}
			}
			  // define the maximal change tolerance based on the aromaticity
			  // of the atom, since we might see "unbalanced" changes due to that
			int maxChange = (isAromaticNode) ? 1 : 0 ;

			  // check all connected edges
			for( boost::tie(aIt, aItEnd)
				= boost::out_edges(*vIt, coreGraph); aIt!=aItEnd; ++aIt)
			{
				  // update electron decrease and increase sum accordingly
				switch (eLevel[*aIt]) {
				case RULE_CONTEXT :
					  // bond unchanged
					break;
				case RULE_LEFT_SIDE :
					  // bond removed
					dec += MoleculeUtil::getBondData( edgeLabel[*aIt] )->valence;
					break;
				case RULE_RIGHT_SIDE :
					  // bond inserted
					inc += MoleculeUtil::getBondData( edgeLabel[*aIt] )->valence;
					break;
				default :
					assert( false /* unknown edge rule context */ );
				}
			}
			  // check if overall electron movement is balanced
			if ( abs( (inc - dec) - chargeChange + protonChange) > maxChange ) {

				std::cerr <<"## DEBUG : checkElectronChange( " <<*vIt <<", " <<nLeftLabel[*vIt] <<(nLevel[*vIt]==RULE_LABEL_CHANGE?", ":"") <<(nLevel[*vIt]==RULE_LABEL_CHANGE?nRightLabel[*vIt]:"") <<" ) = "
							<<" bondLeft="<<dec <<" bondRight="<<inc <<" chargeChange="<<chargeChange <<" protonChange="<<protonChange <<"\n"
						;

				return false;
			}
		}
		  // all ok ...
		return true;
	}

//##############################################################################

	bool
	ChemRule::
	insertGroups( ggl::chem::ChemRule::CoreGraph & core, const ggl::chem::GroupMap & groups ) throw(std::runtime_error)
	{

		using namespace boost;
		using namespace ggl::chem;

		  // access to the property maps to fill
		boost::property_map< CoreGraph, NodeContextProperty >::type
			nodeLevel = boost::get( NodeContextProperty(), core );
		boost::property_map< CoreGraph, NodeLabelProperty >::type
			nodeLabel = boost::get( NodeLabelProperty(), core );
		boost::property_map< CoreGraph, NodeRightLabelProperty >::type
			nodeRightLabel = boost::get( NodeRightLabelProperty(), core );
		boost::property_map< CoreGraph, EdgeLabelProperty >::type
			edgeLabel = boost::get( EdgeLabelProperty(), core );
		boost::property_map< CoreGraph, EdgeContextProperty >::type
			edgeLevel = boost::get( EdgeContextProperty(), core );


		CoreGraph::vertex_descriptor atom, newAtom;
		CoreGraph::edge_descriptor newEdge;

		bool allGroupsReplaced = true;

		  // iterate over all original core nodes (excluding new components)
		  // (avoids recursive group replacement)
		const size_t oldMolSize = boost::num_vertices(core);
		for ( size_t atomID = 0; atomID < oldMolSize; ++atomID ) {

			  // access node
			atom = boost::vertex( atomID, core );

			std::string &atomLabel = nodeLabel[atom];

			if (!ggl::chem::MoleculeUtil::isGroupLabel(atomLabel)) {
				continue;
			}

			  // derive groupID
			std::string groupID = MoleculeUtil::getGroupLabel(atomLabel);
			  // check if groupID known
			if (groups.find(groupID) == groups.end()) {
				allGroupsReplaced = false;
				 // ignore this atom label
				continue;
//				std::ostringstream oss;
//				oss	<<"ggl::chem::ChemRule::insertGroup : the node " <<atomID
//					<<" shows the unknown group label '"<<groupID
//					<<"' within its node label '"<<atomLabel<<"'!";
//				throw std::runtime_error(oss.str());
			}

			  // check if the group node is a context node, only here allowed
			if (nodeLevel[atom] == RULE_LEFT_SIDE || nodeLevel[atom] == RULE_RIGHT_SIDE) {
				std::ostringstream oss;
				oss	<<"ggl::chem::ChemRule::insertGroup : the node " <<atomID
					<<" shows the group label '"<<groupID
					<<"' but is not within rule context or both left-right to encode a charge/class change!";
				throw std::runtime_error(oss.str());
			}

			  // access according group (presence ensured with test from above)
			const MoleculeComponent &group = groups.find(groupID)->second;
			const size_t groupProxy = *group.compIDs.begin();
			  // group graph access
			boost::property_map< ggl::chem::MoleculeComponent::PatternGraph, ggl::PropNodeLabel >::const_type groupNodeLabel = get( ggl::PropNodeLabel(), group.pattern);
			boost::property_map< ggl::chem::MoleculeComponent::PatternGraph, ggl::PropNodeIndex >::const_type groupNodeIndex = get( ggl::PropNodeIndex(), group.pattern);
			boost::property_map< ggl::chem::MoleculeComponent::PatternGraph, ggl::PropEdgeLabel >::const_type groupEdgeLabel = get( ggl::PropEdgeLabel(), group.pattern);


			// relabel compID node (preserve charge/protons = rest of original atom label)
			nodeLabel[atom] = groupNodeLabel[boost::vertex(groupProxy, group.pattern)]
							 + atomLabel.substr(groupID.size(),std::string::npos);

			// check if node is label changing, ie. the right label has to be set too
			if (nodeLevel[atom] == RULE_LABEL_CHANGE) {
				  // get right label
				std::string rightLabel = nodeRightLabel[atom];
				  // ensure right label shows groupID as well
				if (rightLabel.size() < groupID.size() || rightLabel.find(groupID) != 0) {
					std::ostringstream oss;
					oss	<<"ggl::chem::ChemRule::insertGroup : the node " <<atomID
						<<" shows the group label '"<<groupID
						<<"' in left but not within right specification! No such label change possible!";
					throw std::runtime_error(oss.str());
				}
				  // replace group ID with proxy node label and rest of the original label
				nodeRightLabel[atom] = groupNodeLabel[boost::vertex(groupProxy, group.pattern)]
								 + rightLabel.substr(groupID.size(),std::string::npos);
			}

			// insert remaining atoms/bonds of molecule component

			  // mapping of component atom indices to inserted atoms within core
			std::map< MoleculeComponent::PatternGraph::vertex_descriptor
					, CoreGraph::vertex_descriptor> group2core;
			  // proxy node exists already -> mapping known
			group2core[boost::vertex(*group.compIDs.begin(), group.pattern)] = atom;

			  // add and set all group atoms (excluding proxy since already there)
			MoleculeComponent::PatternGraph::vertex_iterator  vi, vi_end;
			for(tie(vi, vi_end)=boost::vertices(group.pattern); vi!=vi_end; ++vi) {
				if (groupNodeIndex[*vi] != groupProxy) {
					  // add new node
					newAtom = boost::add_vertex( core );
					  // store mapping
					group2core[*vi] = newAtom;
					  // set label
					nodeLabel[newAtom] = groupNodeLabel[*vi];
					  // set level to context
					nodeLevel[newAtom] = RULE_CONTEXT;
				}
			}
			  // add and set all group bonds
			MoleculeComponent::PatternGraph::edge_iterator  ei, ei_end;
			for(tie(ei, ei_end)=boost::edges(group.pattern); ei!=ei_end; ++ei) {
				  // add new edge
				newEdge = boost::add_edge( group2core[boost::source(*ei,group.pattern)]
										 , group2core[boost::target(*ei,group.pattern)]
										 , core).first;
				  // set label of new edge
				edgeLabel[newEdge] = groupEdgeLabel[*ei];
				  // set level to context
				edgeLevel[newEdge] = RULE_CONTEXT;
			}
		}

		 // return whether or not all groups have been replaced
		return allGroupsReplaced;

	}

////////////////////////////////////////////////////////////////////////////////

  } // namespace chem
} // namespace ggl
