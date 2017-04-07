
#include "sgm/RP_Hanser96.hh"

#include "ggl/chem/AP_NSPDK.hh"

#include <limits>
#include <cassert>
#include <cmath>

namespace ggl {
namespace chem {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	void
	AP_NSPDK::
	findAllRings( const Molecule & mol )
	{

		  // store all rings up to the maximally predictable size
		sgm::RP_Hanser96 ringFinder;
		  // setup graph interface
		Molecule_Graph molGraph(mol);
		  // run ring perception
		ringFinder.findRings( molGraph, *this, aromaticityModel->maxRingSize );

	}

////////////////////////////////////////////////////////////////////////////////

	void
	AP_NSPDK::
	identifyAromaticEdges( const Molecule & mol )
	{
		assert(aromaticityModel != NULL /*no model selected*/);


		  // setup NSPDK graph instance
		initializeGraph( mol
					, aromaticityModel->ringEdgeLabel
					, aromaticityModel->ringNodeLabelPrefix );

		//////////////////// PREDICT ALL RINGS /////////////////////

		  // handle each ring independently
		for ( size_t curRing = 0; curRing < allRings.size(); ++curRing ) {

			  // test each remaining ring for aromaticity, set the according
			  // flag, and store the aromatic ring edges in aromaticEdges

			const double predValue = getPredictionValue( *(allRings.at(curRing)) );
			allRings[curRing]->predCertainty = (double)std::abs(predValue);

			if (isAromaticRing(predValue)) {
				  // set aromatic
				allRings[curRing]->predState = RingDescriptor::Aromatic;
				  // store ring edges as aromatic ring
				aromaticEdges.insert( allRings.at(curRing)->edges.begin(), allRings.at(curRing)->edges.end() );
			} else {
				  // set non-aromatic
				allRings[curRing]->predState = RingDescriptor::NonAromatic;
			}
		}


	}

////////////////////////////////////////////////////////////////////////////////



	std::vector< std::pair< sgm::RingReporter::RingList, nspdk::SVector> >
	AP_NSPDK::
	getFeatures( const Molecule & mol )
	{

		assert(aromaticityModel != NULL /*no model selected*/);


		  // clear temporary data
		clearData();
		  // identify all rings
		findAllRings( mol );

		  // prune rings to remove fused representation
		AromaticityPerception::pruneFusedRings( allRings );

		  // prune rings that cannot be aromatic
		AromaticityPerception::pruneNonSingleDoubleBondRings( allRings, mol );

		  // setup NSPDK graph instance
		initializeGraph( mol
						, aromaticityModel->ringEdgeLabel
						, aromaticityModel->ringNodeLabelPrefix );

		  // setup container for ring features
		std::vector< std::pair< RingList, nspdk::SVector> > rings(allRings.size());

		  // fill ring data
		for ( size_t curRing = 0; curRing < allRings.size(); ++curRing ) {

			  // store ring list representation
			rings[curRing].first = allRings.at(curRing)->ring;
			  // generate features with model
			rings[curRing].second = getFeaturesOfRing( nspdkGraph, *(allRings.at(curRing)), *aromaticityModel );
		}

		return rings;
	}



////////////////////////////////////////////////////////////////////////////////



	nspdk::SVector
	AP_NSPDK::
	getFeaturesOfRing( nspdk::GraphClass &molGraph
						, const RingDescriptor & ring
						, const AP_NSPDK_Model& model  )
	{
		  // get view points
		sgm::NSPDK_port::ViewPointList ringPoints( ring.ring.begin(), ring.ring.end() );

		  // relabel view point nodes to highlight current ring
		sgm::NSPDK_port::ViewPointList::const_iterator i;
		const std::string& ringViewLabelPrefix = model.ringViewLabelPrefix;
		i=ringPoints.begin(); // ignore first to avoid double labeling
		for (i++; i!=ringPoints.end(); ++i ) {
			  // add ring flag in the beginning
			molGraph.SetVertexLabel( *i, ringViewLabelPrefix + molGraph.GetVertexLabel(*i) );
		}
		  // relabel view point edges to highlight current ring
		sgm::NSPDK_port::ViewPointList::const_iterator l = ringPoints.begin();
		i = ringPoints.begin();
		for (++i; i<ringPoints.end(); ++i, ++l) {
			molGraph.SetEdgeLabel( *l, *i, ringViewLabelPrefix + molGraph.GetEdgeLabel( *l, *i));
			molGraph.SetEdgeLabel( *i, *l, ringViewLabelPrefix + molGraph.GetEdgeLabel( *i, *l));
		}
//std::cerr <<"next ring = ";
//for (RingList::const_iterator x=ring.ring.begin(); x!=ring.ring.end(); ++x) {
//	std::cerr <<*x <<" ";
//}

		 // get list of nodes for which features are to be generated
		sgm::NSPDK_port::ViewPointList viewPoints = getViewPoints(molGraph, ring, model);

		  // get features according to the parameterization of the model that is
		  // used for the prediction afterwards
		nspdk::SVector curFeatures
				= sgm::NSPDK_port::getFeatures(	molGraph
												, model.nspdk_maxDistance
												, model.nspdk_maxRadius
												, viewPoints
												, model.nspdk_featureBitSize
											);

		  // undo view point node relabeling
		const size_t viewLabelLength = ringViewLabelPrefix.size();
		i=ringPoints.begin(); // ignore first to avoid double labeling
		for (i++; i!=ringPoints.end(); ++i ) {
			  // strip view point flag from the beginning
			molGraph.SetVertexLabel( *i, molGraph.GetVertexLabel(*i).substr(viewLabelLength,std::string::npos) );
		}
		  // undo view point edge relabeling
		l = ringPoints.begin();
		i = ringPoints.begin();
		for (++i; i<ringPoints.end(); ++i, ++l) {
			  // strip view point flag from the beginning
			molGraph.SetEdgeLabel( *l, *i, molGraph.GetEdgeLabel( *l, *i).substr(viewLabelLength,std::string::npos) );
			molGraph.SetEdgeLabel( *i, *l, molGraph.GetEdgeLabel( *i, *l).substr(viewLabelLength,std::string::npos) );
		}

		return curFeatures;
	}

////////////////////////////////////////////////////////////////////////////////



	std::vector< std::pair< AP_NSPDK::RingDescriptor, AP_NSPDK::RingWeightVec > >
	AP_NSPDK::
	getNodeWeights( const Molecule & mol )
	{

		assert(aromaticityModel != NULL /*no model selected*/);


		  // clear temporary data
		clearData();
		  // identify all rings
		findAllRings( mol );

		  // prune rings to remove fused representation
		AromaticityPerception::pruneFusedRings( allRings );

		  // prune rings that cannot be aromatic
		AromaticityPerception::pruneNonSingleDoubleBondRings( allRings, mol );

		  // setup NSPDK graph instance
		initializeGraph( mol
						, aromaticityModel->ringEdgeLabel
						, aromaticityModel->ringNodeLabelPrefix );

		  // setup container for ring features
		std::vector< std::pair< RingDescriptor, RingWeightVec > > rings(allRings.size());

		  // fill ring data
		for ( size_t curRing = 0; curRing < allRings.size(); ++curRing ) {

			  // copy ring information for return value
			rings[curRing].first = *(allRings.at(curRing));
			  // update ring type information
			double predictionValue = getPredictionValue( *(allRings.at(curRing)) );
			rings[curRing].first.predState = isAromaticRing(predictionValue) ? RingDescriptor::Aromatic : RingDescriptor::NonAromatic;
			  // set certainty to absolute predictionvalue
			rings[curRing].first.predCertainty = (double)std::abs( predictionValue );

			  // get and store weight list for ring indices
			rings[curRing].second = getWeightsOfRing( nspdkGraph, *(allRings.at(curRing)), *aromaticityModel );
		}

		return rings;
	}

////////////////////////////////////////////////////////////////////////////////



	AP_NSPDK::
	RingWeightVec
	AP_NSPDK::
	getWeightsOfRing( nspdk::GraphClass &molGraph
						, const RingDescriptor & ring
						, const AP_NSPDK_Model& model )
	{
		  // get view points
		sgm::NSPDK_port::ViewPointList ringPoints( ring.ring.begin(), ring.ring.end() );
		assert(ring.ring.size() > 2);

		  // relabel view point nodes to highlight current ring
		sgm::NSPDK_port::ViewPointList::const_iterator i;
		const std::string& ringViewLabelPrefix = model.ringViewLabelPrefix;
		i=ringPoints.begin(); // ignore first to avoid double labeling
		for (i++; i!=ringPoints.end(); ++i ) {
			  // add ring flag in the beginning
			molGraph.SetVertexLabel( *i, ringViewLabelPrefix + molGraph.GetVertexLabel(*i) );
		}
		  // relabel view point edges to highlight current ring
		sgm::NSPDK_port::ViewPointList::const_iterator l = ringPoints.begin();
		i = ringPoints.begin();
		for (++i; i<ringPoints.end(); ++i, ++l) {
			molGraph.SetEdgeLabel( *l, *i, ringViewLabelPrefix + molGraph.GetEdgeLabel( *l, *i));
			molGraph.SetEdgeLabel( *i, *l, ringViewLabelPrefix + molGraph.GetEdgeLabel( *i, *l));
		}
//std::cerr <<"next ring = ";
//for (RingList::const_iterator x=ring.ring.begin(); x!=ring.ring.end(); ++x) {
//	std::cerr <<*x <<" ";
//}

		 // get list of nodes for which features are to be generated
		sgm::NSPDK_port::ViewPointList viewPoints = getViewPoints(molGraph,ring,model);

		  // get features according to the parameterization of the model that is
		  // used for the prediction afterwards
		std::vector<nspdk::SVector> featuresPerViewpoint
				= sgm::NSPDK_port::getNodeFeatures(	molGraph
												, model.nspdk_maxDistance
												, model.nspdk_maxRadius
												, viewPoints
												, model.nspdk_featureBitSize
											);
		assert( featuresPerViewpoint.size() == viewPoints.size() );

		  // undo view point node relabeling
		const size_t viewLabelLength = ringViewLabelPrefix.size();
		for (i = ringPoints.begin()+1; i!=ringPoints.end(); ++i ) {// ignore first to avoid double labeling
			  // strip view point flag from the beginning
			molGraph.SetVertexLabel( *i, molGraph.GetVertexLabel(*i).substr(viewLabelLength,std::string::npos) );
		}
		  // undo view point edge relabeling
		l = ringPoints.begin();
		i = ringPoints.begin();
		for (++i; i<ringPoints.end(); ++i, ++l) {
			  // strip view point flag from the beginning
			molGraph.SetEdgeLabel( *l, *i, molGraph.GetEdgeLabel( *l, *i).substr(viewLabelLength,std::string::npos) );
			molGraph.SetEdgeLabel( *i, *l, molGraph.GetEdgeLabel( *i, *l).substr(viewLabelLength,std::string::npos) );
		}

		  // compute node weights
		RingWeightVec nodeWeights(featuresPerViewpoint.size());
		for (size_t i=0; i<nodeWeights.size(); ++i) {
			  // save ID
			nodeWeights[i].first = viewPoints.at(i);
			  // get the weight for each node
			double curWeight = model.predictValue( featuresPerViewpoint.at(i) );
			  // normalize and store weight
			nodeWeights[i].second = (1-std::exp(-curWeight))/(1+std::exp(-curWeight));
		}

		  // return final node weights
		return nodeWeights;
	}

////////////////////////////////////////////////////////////////////////////////

	double
	AP_NSPDK::
	getPredictionValue( const AromaticityPerception::RingDescriptor & curRing )
	{

		assert(aromaticityModel != NULL /*no model selected*/);

		//////////////////// PREDICT RING /////////////////////

		return aromaticityModel->predictValue(
					getFeaturesOfRing( nspdkGraph, curRing, *aromaticityModel ) );
	}

////////////////////////////////////////////////////////////////////////////////



	bool
	AP_NSPDK::
	isAromaticRing( const AromaticityPerception::RingDescriptor & curRing )
	{

		assert(aromaticityModel != NULL /*no model selected*/);

		//////////////////// PREDICT RING /////////////////////

		return isAromaticRing( getPredictionValue( curRing ));
	}

////////////////////////////////////////////////////////////////////////////////



	AP_NSPDK::
	ViewPointList
	AP_NSPDK::
	getViewPoints( nspdk::GraphClass &molGraph
					, const RingDescriptor & ring
					, const AP_NSPDK_Model & model ) const
	{
		ViewPointList viewPoints;

		if (model.ringCentered) {

			if (model.nspdk_maxRingDistance < 1) {
				  // add all but the first node since it is also the last node
				RingList::const_iterator ri = ring.ring.begin();
				if (ri != ring.ring.end()) { ++ri; }
				viewPoints = ViewPointList( ri, ring.ring.end() );
			} else {
				  // init list of nodes to be checked in next layer
				std::set< size_t > toCheck(ring.ring.begin(),ring.ring.end());
				  // this is the list of all known nodes in the last layers
				std::set< size_t > lastLayers(ring.ring.begin(),ring.ring.end());
				  // this will be filled for each next layer
				std::set< size_t > nextLayer;
				  // iteratively generate layer after layer starting from the ring nodes
				for (size_t layer = 1; layer <= model.nspdk_maxRingDistance; ++layer) {
					  // for each node to be checked
					while (!toCheck.empty()) {
						  // pop first element
						size_t curNode = *(toCheck.begin());
						toCheck.erase(toCheck.begin());
						// generate neighborhood of curNode
						std::vector<unsigned> neighbors = molGraph.GetFixedDistanceVertexIDList( (unsigned)curNode, 1 );
						for (std::vector<unsigned>::const_iterator neigh = neighbors.begin(); neigh != neighbors.end(); ++neigh)  {
							// check if neighbor is unknown
							if (lastLayers.find(*neigh) == lastLayers.end()
									&& nextLayer.find(*neigh) == nextLayer.end() )
							{
								// add to nextLayer
								nextLayer.insert(*neigh);
							}
						}
					}
					  // make next layer to be checked
					toCheck = nextLayer;
					  // mark as known
					lastLayers.insert(toCheck.begin(), toCheck.end());
					  // clear data for next iteration
					nextLayer.clear();
				}
				  // generate final view point list
				viewPoints = ViewPointList( lastLayers.begin(), lastLayers.end() );
			}
		} else {
			  // add all nodes
			viewPoints = ViewPointList( molGraph.VertexSize(), 0 );
			for (size_t x=0; x<viewPoints.size(); ++x) {
				viewPoints[x] = x;
			}
		}

		return viewPoints;
	}
////////////////////////////////////////////////////////////////////////////////



	bool
	AP_NSPDK::
	isAromaticRing( const double predValue )
	{

		assert(aromaticityModel != NULL /*no model selected*/);

		//////////////////// PREDICT RING /////////////////////

		return aromaticityModel->predict( predValue );
	}

////////////////////////////////////////////////////////////////////////////////


	void
	AP_NSPDK::
	initializeGraph( const Molecule mol
					, const std::string edgeLabel
					, const std::string ringFlag )
	{

		  // create a dummy graph to create NSPDK graph instance
		Molecule_Graph molGraph(mol);

		  // reset NSPDK graph instance for feature generation
		nspdkGraph = nspdk::GraphClass();
		  // initialize NSPDK graph with the molecule graph
		  // (still relabeling of rings and protons needed)
		sgm::NSPDK_port::setGraphFromGraphInterface( molGraph, nspdkGraph );


		  // temporary variables
		const size_t ringFlagLength = ringFlag.size();
		const MoleculeUtil::BondLabelData * bondData = NULL;

		  // relabel all rings within nspdkGraph
		for (size_t curRing=0; curRing<allRings.size(); ++curRing) {

			const RingList & ringList = allRings.at(curRing)->ring;

			  // check if ring contains only single or double bonds
			RingList::const_iterator l = ringList.begin();
			RingList::const_iterator i = ringList.begin();
#ifndef NDEBUG
			for ( ++i; i != ringList.end(); ++i,++l) {
				const std::string& label = nspdkGraph.GetEdgeLabel( *l, *i );
				bondData = MoleculeUtil::getBondData( label );
				if( bondData != NULL ) {
					  // ensure bond is possible aromatic ring part, i.e. valence == (1 or 2)
					assert( bondData->valence <= 2 );
					assert( bondData->valence >= 1 );
				}
			}
#endif
			  // relabel all edges (both directions)
			l = ringList.begin();
			i = ringList.begin();
			for ( ++i; i != ringList.end(); ++i,++l) {
				nspdkGraph.SetEdgeLabel( *l, *i, edgeLabel );
				nspdkGraph.SetEdgeLabel( *i, *l, edgeLabel );
			}
			  // relabel all nodes
			for ( i=ringList.begin(); i!=ringList.end(); ++i ) {
				  // check if node label already set

				std::string label = nspdkGraph.GetVertexLabel( *i );
				if ( label.size() <= ringFlagLength || ringFlag != label.substr(0,ringFlagLength) ) {
					  // check if currently aromatic : if so get non-aromatic label
					if (MoleculeUtil::getAtomData(MoleculeUtil::getAtom( label ))->isAromatic != 0) {
						label = MoleculeUtil::getComplexAtomLabel(
								*MoleculeUtil::getAromaticPendant( MoleculeUtil::getAtom( label ) )
								, MoleculeUtil::getProtons(label)
								, MoleculeUtil::getCharge(label)
								, MoleculeUtil::getClass(label)
								, false
								);
					}
					  // set "uncertain" label
					nspdkGraph.SetVertexLabel( *i, ringFlag+label );
				}
			}

		} // all rings
	}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}} // namespace
