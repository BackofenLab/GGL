/* -*- mode:c++ -*- */
#ifndef NSPDK_NSPDK_FEATUREGENERATOR_H
#define NSPDK_NSPDK_FEATUREGENERATOR_H

#include "nspdk/Utility.h"
#include "nspdk/FlagsService.h"
#include "nspdk/GraphClass.h"
#include "nspdk/FeatureGenerator.h"
#include <limits>
//#include <omp.h>

using namespace std;

namespace nspdk {

class DebugClass {
public:
	map<unsigned,string> mHashToFeatureMap;
	map<unsigned,string> mHashToGraphMap;
public:
	void Clear();
	void StoreFeatureCodeToFeatureInfo(unsigned aFeatureCode, vector<unsigned>& aDetailsList);
	void SerializedRootedGraphCanonicalFormEncoding(unsigned aDiscreteEncoding, int aRootVertexIndex, const GraphClass& aG, int aRadius);
	void Output(ostream& out) const;
	void OutputFeatureEncoding(ostream& out) const;
};

//----------------------------------------------------------------------------------------------------------------------------------
class NSPDK_FeatureGenerator: public FeatureGenerator, public FlagsServiceClient {
protected:
	unsigned mRadius;
	unsigned mDistance;
	string mMatchType;
	unsigned mHashBitSize;
	unsigned mHashBitMask;
	bool mMinKernel;
	bool mNormalization;
	bool mUseRealVectorInformation;
	unsigned mNumHashFunctionsMinHash;
	unsigned mVertexDegreeThreshold;
	unsigned mDebugVerbosity;
	vector<vector<unsigned> > mFeatureCache;
	mutable DebugClass mDebugInfo;

public:
	NSPDK_FeatureGenerator(const std::string& id);
	virtual ~NSPDK_FeatureGenerator(){}
	void OutputParameters(ostream& out) const;
	void OutputFeatureMap(ostream& out) const;
	virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	virtual void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	unsigned GenerateVectorHashCode(SVector& x);
	unsigned GenerateVertexHashCode(unsigned aSrcID, const GraphClass& aG);
	virtual void GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);
	virtual void GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, SVector& x);
	void InitFeatureCache(const GraphClass& aG, unsigned aRadius);
	unsigned GenerateVertexNeighbourhoodHashCode(unsigned aRootID, const GraphClass& aG, unsigned aRadius);
	unsigned GenerateGraphHashCode(const GraphClass& aG, const vector<unsigned>& aFirstEndpointList);
	void ConvertSparseVectorToMinFeatureVector(SVector& x);
	void CanonicalGraphVertexList(const GraphClass& aG, vector<unsigned>& oVertexList);
	unsigned GraphCanonicalFormEncoding(const GraphClass& aG);
	unsigned RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius);
	unsigned MultiRootedGraphCanonicalFormEncoding(vector<int>& aRootVertexIndexList, const GraphClass& aG);
	unsigned Radius0RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG);
	unsigned Radius1RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG);
	unsigned RadiusKRootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius);
	unsigned RootedGraphCanonicalFormEncoding(const GraphClass& aG, unsigned aRootID);
	void GetFirstEndpoints(const GraphClass& aG, vector<unsigned>& oFirstEndpointRootId) const;
	virtual unsigned HashFunc(const string& aString, unsigned aBitMask = 2147483647);
	virtual unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask = 2147483647);
};

//----------------------------------------------------------------------------------------------------------------------------------
//Shell Kernel
class SK_FeatureGenerator: public NSPDK_FeatureGenerator {
protected:
	vector<vector<vector<unsigned> > > mShellCache;
public:
	SK_FeatureGenerator(const std::string& id);
	virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	virtual void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	void GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);
	unsigned OuterRadiusKInnerRadiusMGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aOuterRadius, int aInnerRadius);
	void InitShellCache(const GraphClass& aG, unsigned aRadius);
	unsigned GenerateOuterRadiusKInnerRadiusMVertexNeighbourhoodHashCode(int aRootVertexIndex, const GraphClass& aG, int aOuterRadius, int aInnerRadius);
};

//----------------------------------------------------------------------------------------------------------------------------------
//Union of the Shortest Paths Kernel
class USPK_FeatureGenerator: public NSPDK_FeatureGenerator {
public:
	USPK_FeatureGenerator(const std::string& id);
	void GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);
};

//----------------------------------------------------------------------------------------------------------------------------------
//Weighted Decomposition Kernel
class WDK_FeatureGenerator: public NSPDK_FeatureGenerator {
public:
	WDK_FeatureGenerator(const std::string& id);
	virtual ~WDK_FeatureGenerator() {
	}
	void InitFeatureCache(const GraphClass& aG);
	virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList);
	virtual void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	virtual void GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);
	void GenerateVertexFeatures(unsigned aRootVertexIndex, const GraphClass& aG, unsigned aRadius, SVector& x);
	void ReHash(SVector& x, unsigned aReHashCode);
protected:
	unsigned mLowerVertexDegreeThreshold;
	vector<SVector> mSparseFeatureCache;
	vector<bool> mSparseFeatureFlagCache;
};

//----------------------------------------------------------------------------------------------------------------------------------
//Abstract NSPDK
class ANSPDK_FeatureGenerator: public NSPDK_FeatureGenerator {
public:
	ANSPDK_FeatureGenerator(const std::string& id);
	virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	virtual void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
	void GenerateAbstractVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);
};

//----------------------------------------------------------------------------------------------------------------------------------
//Poly Bias Kernel
class PBK_FeatureGenerator: public NSPDK_FeatureGenerator {
public:
	PBK_FeatureGenerator(const std::string& id);
	virtual ~PBK_FeatureGenerator() {
	}
	virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList);
	virtual void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
public:
	ANSPDK_FeatureGenerator mANSPDK;
	WDK_FeatureGenerator mWDK;
};

} // namespace nsdpk

#endif /* NSPDK_NSPDK_FEATUREGENERATOR_H */

