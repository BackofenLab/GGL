/*
 * NSPDK_port.hh
 *
 *  Created on: 11.05.2011
 *      Author: mmann
 */

#ifndef NSPDK_PORT_HH_
#define NSPDK_PORT_HH_

#include "sgm/Graph_Interface.hh"
#include "nspdk/GraphClass.h"
#include "nspdk/NSPDK_FeatureGenerator.h"


namespace sgm {

	 /*! @brief Ports NSPDK functionality to SGM
	  *
	  * A utility class that ports the functionalities of the NSPDK graph kernel
	  * for its application the SGM library.
	  *
	  * @author Martin Mann (c) 2008 http://www.bioinf.uni-freiburg.de/~mmann/
	  */
	class NSPDK_port {

	public:

		  //! Sparse vector definition that represents the features of a graph
		typedef nspdk::SVector FeatureVector;

		  //! List of view points to localize the feature generation
		typedef std::vector< unsigned > ViewPointList;

	protected:

		 //! The standard NSPDK feature generator
		static
		nspdk::NSPDK_FeatureGenerator featureGenerator_nspdk;

	public:

		 /*!
		  * Initializes an NSPDK graph object from a given Graph_Interface.
		  *
		  * @param iG IN the graph to be used for the initialization
		  * @param oG OUT the EMPTY graph object to be filled
		  *
		  */
		static
		void setGraphFromGraphInterface(	const Graph_Interface & iG
											, nspdk::GraphClass & oG);

		 /*!
		  * Computes the feature vector for a given graph
		  *
		  * @param g IN the graph that has to be screened for features
		  * @param maxDistance the maximal distance between subgraphs to
		  *        generate features
		  * @param maxRadius the maximal radius between subgraphs to
		  *        generate features
		  * @param viewPoints IN optional list of view points to localize the
		  *        feature generation to the according subgraph and its
		  *        surrounding
		  * @param featureBitSize number of bits to represent the feature space
		  *
		  * @return the feature vector describing the graph g
		  */
		static
		FeatureVector
		getFeatures(	nspdk::GraphClass & g
						, const size_t maxDistance
						, const size_t maxRadius
						, const ViewPointList & viewPoints = ViewPointList()
						, const int featureBitSize = 30
						);


		 /*!
		  * Computes the feature vector for each node of a given graph
		  *
		  * @param g IN the graph that has to be screened for features
		  * @param maxDistance the maximal distance between subgraphs to
		  *        generate features
		  * @param maxRadius the maximal radius between subgraphs to
		  *        generate features
		  * @param viewPoints IN optional list of view points to localize the
		  *        feature generation to the according subgraph and its
		  *        surrounding
		  * @param featureBitSize number of bits to represent the feature space
		  *
		  * @return the feature vector describing the graph g for each node
		  */
		static
		std::vector<NSPDK_port::FeatureVector>
		getNodeFeatures(	nspdk::GraphClass & g
							, const size_t maxDistance
							, const size_t maxRadius
							, const ViewPointList & viewPoints = ViewPointList()
							, const int featureBitSize = 30
						);


	};


} // namespace sgm



#endif /* NSPDK_PORT_HH_ */
