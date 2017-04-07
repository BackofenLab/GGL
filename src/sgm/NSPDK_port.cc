
#include "sgm/NSPDK_port.hh"
#include <cassert>

namespace sgm {

	//////////////////////////////////////////////////////////////////////////

	nspdk::NSPDK_FeatureGenerator
	NSPDK_port::
	featureGenerator_nspdk
		= nspdk::NSPDK_FeatureGenerator(std::string("nspdk"));


	//////////////////////////////////////////////////////////////////////////

	void
	NSPDK_port::
	setGraphFromGraphInterface(	const Graph_Interface &iG, nspdk::GraphClass& oG)
	{
		assert(oG.IsEmpty() /*we assume an empty graph*/);

		using namespace nspdk;
		using namespace std;

		  // prepare node information to store
		vector<bool> vertex_status;
			vertex_status.push_back(true);//kernel point
			vertex_status.push_back(true);//kind
			vertex_status.push_back(true);//viewpoint
			vertex_status.push_back(false);//dead
			vertex_status.push_back(false);//forbidden distance
			vertex_status.push_back(false);//forbidden neighborhood
		  // will hold the label information (updated for each node)
		vector<string> vertex_symbolic_attribute_list(1,string("NO_LABEL_SET"));

		// copy all nodes
		for (size_t i=0; i<iG.getNodeNumber(); ++i) {
			  // create node
			oG.InsertVertex();
			  // update and store label information
			vertex_symbolic_attribute_list[0] = iG.getNodeLabel(i);
			oG.SetVertexSymbolicAttributeList( i, vertex_symbolic_attribute_list);
			  // set status information (equal for all new nodes)
			oG.SetVertexStatusAttributeList( i, vertex_status );
		}

		  // prepare edge information to store
		  // will hold the label information (updated for each edge)
		vector<string> edge_symbolic_attribute_list(1,string("NO_LABEL_SET"));
//		  // will hold addition information (none needed so far)
//		vector<double> edge_numeric_attribute_list(1,0.0);

		// copy all edges
		for (size_t i=0; i < iG.getNodeNumber(); ++i) {
			for (	Graph_Interface::OutEdge_iterator e = iG.getOutEdgesBegin(i),
					eEnd = iG.getOutEdgesEnd(i); e != eEnd; ++e )
			{
				  // insert edge
				unsigned edge_index=oG.InsertEdge(e->getFromIndex(),e->getToIndex());
				  // update and store edge label information
				edge_symbolic_attribute_list[0] = e->getEdgeLabel();
				oG.SetEdgeSymbolicAttributeList( edge_index, edge_symbolic_attribute_list);
//				  // store additional information (none needed so far)
//				oG.SetEdgeNumericAttributeList( edge_index, edge_numeric_attribute_list);

			}
		}

	}

	//////////////////////////////////////////////////////////////////////////

	NSPDK_port::
	FeatureVector
	NSPDK_port::
	getFeatures(	nspdk::GraphClass & g
					, const size_t maxDistance
					, const size_t maxRadius
					, const ViewPointList & viewPoints
					, const int featureBitSize
				)
	{

		  // set according maximal distance and radius
		featureGenerator_nspdk.set_flag( "distance", nspdk::stream_cast<std::string>(maxDistance));
		featureGenerator_nspdk.set_flag( "radius", nspdk::stream_cast<std::string>(maxRadius));

		  // set standard flags
		featureGenerator_nspdk.set_flag("match_type","hard");
		featureGenerator_nspdk.set_flag("hash_bit_size",nspdk::stream_cast<std::string>(featureBitSize));
		featureGenerator_nspdk.set_flag("hash_bit_mask",nspdk::stream_cast<std::string>((2 << (featureBitSize - 1)) - 1));
		featureGenerator_nspdk.set_flag("verbosity","0");

		featureGenerator_nspdk.set_flag("vertex_degree_threshold","10");
		featureGenerator_nspdk.set_flag("num_hash_functions_minhash","400");

		  // the feature vector to fill
		FeatureVector features;
		  // generate features via the set feature generator
		featureGenerator_nspdk.generate_feature_vector( g, features, viewPoints );

		return features;
	}

	//////////////////////////////////////////////////////////////////////////

	std::vector<NSPDK_port::FeatureVector>
	NSPDK_port::
	getNodeFeatures(	nspdk::GraphClass & g
						, const size_t maxDistance
						, const size_t maxRadius
						, const ViewPointList & viewPoints
						, const int featureBitSize
					)
	{

		  // set according maximal distance and radius
		featureGenerator_nspdk.set_flag( "distance", nspdk::stream_cast<std::string>(maxDistance));
		featureGenerator_nspdk.set_flag( "radius", nspdk::stream_cast<std::string>(maxRadius));

		  // set standard flags
		featureGenerator_nspdk.set_flag("match_type","hard");
		featureGenerator_nspdk.set_flag("hash_bit_size",nspdk::stream_cast<std::string>(featureBitSize));
		featureGenerator_nspdk.set_flag("hash_bit_mask",nspdk::stream_cast<std::string>((2 << (featureBitSize - 1)) - 1));
		featureGenerator_nspdk.set_flag("verbosity","0");

		featureGenerator_nspdk.set_flag("vertex_degree_threshold","10");
		featureGenerator_nspdk.set_flag("num_hash_functions_minhash","400");


		  // the feature vectors to fill
		std::vector<FeatureVector> features;
		  // generate features via the set feature generator
		if (!viewPoints.empty()) {
			  // generate features just for requested list
			featureGenerator_nspdk.generate_vertex_feature_vector( g, features, viewPoints );
		} else {
			  // generate features just for all nodes
			ViewPointList allPoints(g.VertexSize());
			for (size_t i=0; i<allPoints.size(); ++i) {
				allPoints[i] = i;
			}
			featureGenerator_nspdk.generate_vertex_feature_vector( g, features, allPoints );
		}

		return features;
	}

	//////////////////////////////////////////////////////////////////////////

} // namespace
