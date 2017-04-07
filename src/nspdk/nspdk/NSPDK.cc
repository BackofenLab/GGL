#include <ctime>
#include <numeric>


#include "nspdk/Utility.h"
#include "nspdk/SVector.h"
#include "nspdk/Histogram.h"
#include "nspdk/BaseGraphClass.h"
#include "nspdk/GraphClass.h"
#include "nspdk/NSPDK_FeatureGenerator.h"


using namespace std;
using namespace nspdk;

//---------------------------------------------------------------------------------
const string CREDITS("Name: NeighborhoodSubgraphPaiwiseDistanceKernel\nVersion: 7.3\nProgrammer: Fabrizio Costa\nLast update: March 1 2011");
//---------------------------------------------------------------------------------
ofstream LOGF("log",std::ios_base::app);

//FlagsService& The_FlagsService= FlagsService::get_instance();

inline double random01(){
  return (double)rand()/(double)RAND_MAX;
}

inline int IntHash(int key, int aModulo=RAND_MAX)
{
  key = ~key + (key << 15); // key = (key << 15) - key - 1;
  key = key ^ (key >> 12);
  key = key + (key << 2);
  key = key ^ (key >> 4);
  key = key * 2057; // key = (key + (key << 3)) + (key << 11);
  key = key ^ (key >> 16);
  return key%aModulo;
}

class Timer {
public:
  typedef double diff_type;

  // Same as Timer t; t.begin();
  Timer(): start(std::clock()), elapsed(0) {}
  ~Timer(){cout<<"Elapsed time: "<<end()<<endl;}
  // Last result before a call to begin()
  diff_type last() const { return elapsed; }
  // Reset the timer
  void begin() { start = std::clock(); elapsed = 0; }
  // Save the result
  diff_type end(){
    elapsed = (diff_type)std::clock() - start;
    elapsed /= CLOCKS_PER_SEC;
    return elapsed;
  }

private:
  std::clock_t start;
  diff_type    elapsed;
};

class ProgressBar{
 public:
  ProgressBar():mCounter(0){}
  ~ProgressBar(){cout<<endl<<"Counted "<<mCounter<<" times."<<endl;}
  void Begin(){mCounter=0;}
  void Count(){
    mCounter++;
    if (mCounter%10==0) cout<<"."<<flush;
    if (mCounter%100==0) cout<<"|"<<flush;
    if (mCounter%1000==0) cout<<mCounter/(1000)<<"K"<<flush;
    if (mCounter%1000000==0) cout<<mCounter/(1000000)<<"M"<<flush;
  }
  unsigned End(){return mCounter;}  
 private:
  unsigned mCounter;
  Timer mTimer;
};

class ParameterWrapperClass{
public:
  ParameterWrapperClass():mGspanInputFileName(""),
			  mSparseBinaryInputFileName(""),
			  mSparseASCIIInputFileName(""),
			  mSparseASCIITestInputFileName(""),
			  mTrainTargetInputFileName(""),
			  mRadiusMax(1),
			  mDistanceMax(4),
			  mMatchingType("hard"),
			  mVerbose(false),
			  mFeatureBitSize(30),
			  mMinKernel(false),
			  mType("nspdk"),
			  mNormalization(true),
			  mNumNearestNeighbors(10),
			  mNumHashFunctions(30),
			  mSampleSize(10),
			  mNonRedundantFilter(1),
			  mAbstractFeatures(false),
			  mOutputFeatures(false),
			  mOutputFeatureMap(false),
			  mOutputKernel(false),  
			  mOutputApproximateKNN(false),  
			  mOutputFastKNN(false),  
			  mOutputTrueKNN(false),  
			  mOutputCluster(false),  
			  mOutputApproximateKNNPrediction(false),
			  mOutputTrueKNNPrediction(false),
			  mOutputDatasetSplit(false),
			  mOutputHashEncoding(false),
			  mCaching(true),
			  mTrueSort(true),
			  mEccessNeighbourSizeFactor(0),
			  mHashFactor(1),
			  mNumCenters(10),
			  mSizeThreshold(30),
			  mImbalanceTolerance(1.5),
			  mMixedMatchingType(false),
			  mDebug(0){}
  
  void Usage(string aCommandName){
    cerr<<"Usage: "<<aCommandName<<endl
	<<"-fg <file name gspan format> for input from file "<<endl
	<<"-fsb <file name sparse format binary> for input from file"<<endl
	<<"-fsa <file name sparse format ascii> for input from file"<<endl
	<<"-fsats <file name sparse format ascii> for input from file for test set"<<endl
	<<"-ftrt <file name> for target for train set"<<endl
	<<"[-of flag to output feature encoding (default: "<< mOutputFeatures<<")]"<<endl
	<<"[-ofm flag to output feature map encoding (default: "<< mOutputFeatureMap<<")]"<<endl
	<<"[-ok flag to output kernel matrix (default: "<< mOutputKernel<<")]"<<endl
	<<"[-oaknn flag to output approximate k-nearest neighburs (default: "<< mOutputApproximateKNN<<")]"<<endl
	<<"[-ofknn flag to output fast k-nearest neighburs (default: "<< mOutputFastKNN<<")]"<<endl
	<<"[-otknn flag to output true (i.e. implies full kernel matrix evaluation) k-nearest neighburs (default: "<< mOutputTrueKNN<<")]"<<endl	
	<<"[-oc flag to output clusters (default: "<< mOutputCluster<<")]"<<endl
	<<"[-oaknnp flag to output approximate knn prediction (default: "<< mOutputApproximateKNNPrediction<<")]"<<endl
	<<"[-otknnp flag to output approximate knn prediction (default: "<< mOutputTrueKNNPrediction<<")]"<<endl
	<<"[-ods flag to output both training set and test set split (default: "<< mOutputDatasetSplit<<")]"<<endl
	<<"[-ohe flag to output hash encoding (default: "<< mOutputHashEncoding<<")]"<<endl
	<<"[-b <feature space bits size> (default: "<< mFeatureBitSize <<")]"<<endl
	<<"[-R <max radius> (default: "<< mRadiusMax <<")]"<<endl
	<<"[-D <max distance relations> (default: "<< mDistanceMax <<")]"<<endl
	<<"[-af flag to activate abstract features (default: "<< mAbstractFeatures <<")]"<<endl
	<<"[-T <nspdk|anspdk|mnspdk|nspdk3d|alnspdk|nspdkvp> (default: "<< mType <<")]"<<endl
	<<"[-t <hard | soft | hard_soft> as neighbourhood matching use HARD=exact matching, use SOFT=attribute matching with root identifier as radius 0 neighborhood, use HARD-SOFT=attribute matching with root identifier as full neighbourhood encoding (default: "<< mMatchingType <<")]"<<endl
	<<"[-tm  flag to activate mixed neighbourhood matching (default: "<< mMixedMatchingType <<")]"<<endl
	<<"[-mink flag to set minimum kernel rather than dot product (default: "<< mMinKernel <<")]"<<endl
	<<"[-nn flag to de-acivate normalization (default: "<< mNormalization <<")]"<<endl
	<<"[-nhf <num hash functions> for the Locality Sensitive Hashing function (default: "<< mNumHashFunctions <<")]"<<endl
	<<"[-hf <hash factor> for the Locality Sensitive Hashing function (default: "<< mHashFactor <<")]"<<endl
	<<"[-knn <num nearest neighbors> (default: "<< mNumNearestNeighbors<<")]"<<endl
	<<"[-ensf <eccess neighbour size factor> (default: "<< mEccessNeighbourSizeFactor<<") (0 to avoid trimming)] "<<endl
	<<"[-ss <sample size> for clustering procedure (default: "<< mSampleSize <<")]"<<endl
	<<"[-nrt <similarity filtering of redundant centers [0,1]> for clustering procedure (default: "<< mNonRedundantFilter <<") (the smaller the less similar the centers)]"<<endl
	<<"[-no-cache flag to deactivate caching of kernel value computation (to minimize memory usage) (default: "<< mCaching <<")]"<<endl
	<<"[-no-true-sort flag to deactivate sorting approximate neighbours with true kernel computation (default: "<< mTrueSort <<")]"<<endl
	<<"[-nc num centers (default: "<< mNumCenters <<")]"<<endl
	<<"[-st <size threshold> (default: "<< mSizeThreshold <<")]"<<endl
	<<"[-it <imbalance tolerance> (default: "<< mImbalanceTolerance <<")]"<<endl
	<<"[-debug <debug level> (default: "<< mDebug <<")]"<<endl;

    exit(0);
  }

  void Init(int argc, char** argv){
    vector<string> options;
    for (int i=1;i<argc; ++i ) options.push_back(argv[i]);
    for (vector<string>::iterator it=options.begin();it!=options.end();++it)
      {
	if ((*it)=="-h" || (*it)=="--help") Usage(argv[0]); 
	else if ((*it)=="-fg") mGspanInputFileName=(*(++it)); 
	else if ((*it)=="-fsb") mSparseBinaryInputFileName=(*(++it)); 
	else if ((*it)=="-fsa") mSparseASCIIInputFileName=(*(++it)); 
	else if ((*it)=="-fsats") mSparseASCIITestInputFileName=(*(++it)); 
	else if ((*it)=="-ftrt") mTrainTargetInputFileName=(*(++it)); 
	else if ((*it)=="-knn") mNumNearestNeighbors=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-ensf") mEccessNeighbourSizeFactor=stream_cast<double>(*(++it)); 
	else if ((*it)=="-hf") mHashFactor=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-R") mRadiusMax=stream_cast<double>(*(++it)); 
	else if ((*it)=="-D") mDistanceMax=stream_cast<double>(*(++it)); 
	else if ((*it)=="-t") mMatchingType=(*(++it));
	else if ((*it)=="-mt") mMixedMatchingType=true;
	else if ((*it)=="-v") mVerbose=true;
	else if ((*it)=="-mink") mMinKernel=true;
	else if ((*it)=="-b") mFeatureBitSize=stream_cast<int>(*(++it)); 
	else if ((*it)=="-T") mType=(*(++it)); 
	else if ((*it)=="-of") mOutputFeatures=true;
	else if ((*it)=="-ofm") mOutputFeatureMap=true;
	else if ((*it)=="-ok") mOutputKernel=true;
	else if ((*it)=="-oaknn") mOutputApproximateKNN=true;
	else if ((*it)=="-ofknn") mOutputFastKNN=true;
	else if ((*it)=="-otknn") mOutputTrueKNN=true;
	else if ((*it)=="-oaknnp") mOutputApproximateKNNPrediction=true;
	else if ((*it)=="-otknnp") mOutputTrueKNNPrediction=true;
	else if ((*it)=="-ods") mOutputDatasetSplit=true;
	else if ((*it)=="-ohe") mOutputHashEncoding=true;
	else if ((*it)=="-oc") mOutputCluster=true;
	else if ((*it)=="-nn") mNormalization=false;
	else if ((*it)=="-debug") mDebug=stream_cast<int>(*(++it)); 
	else if ((*it)=="-nhf") mNumHashFunctions=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-ss") mSampleSize=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-nrt") mNonRedundantFilter=stream_cast<double>(*(++it)); 
	else if ((*it)=="-af") mAbstractFeatures=true;
	else if ((*it)=="-no-cache") mCaching=false;
	else if ((*it)=="-no-true-sort") mTrueSort=false;
	else if ((*it)=="-nc") mNumCenters=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-st") mSizeThreshold=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-it") mImbalanceTolerance=stream_cast<double>(*(++it)); 
	else {cerr<<"Unrecognized parameter: "<<(*it)<<"."<<endl;throw exception();}
      }
    if (!(mMatchingType=="hard"||mMatchingType=="soft"||mMatchingType=="hard_soft"||mMatchingType=="multiview"||mMatchingType=="mixed")) {cerr<<"Wrong value for parameter: -t: "<<mMatchingType<<endl;throw exception();}
  }

public:
  string mGspanInputFileName;
  string mSparseBinaryInputFileName;
  string mSparseASCIIInputFileName;
  string mSparseASCIITestInputFileName;
  string mTrainTargetInputFileName;
  double mRadiusMax;
  double mDistanceMax;
  string mMatchingType;
  bool mVerbose;
  int mFeatureBitSize;
  bool mMinKernel;
  string mType;
  bool mNormalization;
  unsigned mNumNearestNeighbors;
  unsigned mNumHashFunctions;
  unsigned mSampleSize;
  double mNonRedundantFilter;
  bool mAbstractFeatures;
  bool mOutputFeatures;
  bool mOutputFeatureMap;
  bool mOutputKernel; 
  bool mOutputApproximateKNN; 
  bool mOutputFastKNN; 
  bool mOutputTrueKNN; 
  bool mOutputCluster; 
  bool mOutputApproximateKNNPrediction;
  bool mOutputTrueKNNPrediction;
  bool mOutputDatasetSplit;
  bool mOutputHashEncoding;
  bool mCaching;
  bool mTrueSort;
  double mEccessNeighbourSizeFactor;
  unsigned mHashFactor;
  unsigned mNumCenters;
  unsigned mSizeThreshold;
  double mImbalanceTolerance;
  bool mMixedMatchingType;
  int mDebug;
};


class MainProcessClass{
protected:
  NSPDK_FeatureGenerator* pmFeatureGenerator;
  NSPDK_FeatureGenerator* pmSoftFeatureGenerator;
  double mAlpha;
  unsigned mNumHashFunctions;
  unsigned mHashFactor;
  unsigned mNeighborhoodSize;
  unsigned mSampleSize;
  double mNonRedundantFilter;
  bool mAbstractFeatures;
  string mMode;
  bool mCaching;
  bool mTrueSort;
  double mEccessNeighbourSizeFactor;
  vector<SVector> mDataset;  
  vector<multimap<unsigned,unsigned> > mBinDataStructure;
  map<pair<unsigned,unsigned> , double> mKernelMap;
  multimap <unsigned,unsigned> mInvertedIndex;
  map<unsigned,vector<unsigned> > mSignatureMap;
  bool mMixedMatchingType;
  bool mOutputHashEncoding;

public:
  MainProcessClass(NSPDK_FeatureGenerator* paFeatureGenerator, 
		   NSPDK_FeatureGenerator* paSoftFeatureGenerator, 
		   unsigned aNumHashFunctions, 
		   unsigned aHashFactor,
		   unsigned aNeighborhoodSize, 
		   unsigned aSampleSize, 
		   double aNonRedundantFilter, 
		   bool aAbstractFeatures, 
		   string aMode, 
		   bool aCaching,
		   bool aTrueSort,
		   double aEccessNeighbourSizeFactor,
		   bool aMixedMatchingType,
		   bool aOutputHashEncoding):
    pmFeatureGenerator(paFeatureGenerator), 
    pmSoftFeatureGenerator(paSoftFeatureGenerator), 
    mNumHashFunctions(aNumHashFunctions),
    mHashFactor(aHashFactor),
    mNeighborhoodSize(aNeighborhoodSize),
    mSampleSize(aSampleSize),
    mNonRedundantFilter(aNonRedundantFilter),
    mAbstractFeatures(aAbstractFeatures),
    mMode(aMode),
    mCaching(aCaching),
    mTrueSort(aTrueSort),
    mEccessNeighbourSizeFactor(aEccessNeighbourSizeFactor),
    mMixedMatchingType(aMixedMatchingType),
    mOutputHashEncoding(aOutputHashEncoding){}

  void Generate(const GraphClass& aG, SVector& oX){
    pmFeatureGenerator->generate_feature_vector(aG,oX);

    if (mAbstractFeatures){    
      GraphClass abstract_g(aG);
      //remove vertex label information to make abstract graph i.e. graph based on pure connectivity information
      for (unsigned i=0;i<abstract_g.VertexSize();++i){
	abstract_g.SetVertexSymbolicAttributeList(i,0,"na");
      }
      SVector abstract_x;
      pmFeatureGenerator->generate_feature_vector(abstract_g,abstract_x);
      //add features to original ones
      oX.add(abstract_x);
      //re-normalize the overall feature representation
      double norm=dot(oX,oX);
      oX.scale(1/norm);
    }

    if (mMixedMatchingType){
      SVector soft_x;
      pmSoftFeatureGenerator->generate_feature_vector(aG,soft_x);
      //add features to original ones
      oX.add(soft_x);
      //re-normalize the overall feature representation
      double norm=dot(oX,oX);
      oX.scale(1/norm);
    }    
  }

  void InputStringList(const string& aFileName, vector<string>& oStringList){
    cout<<"Reading "<<aFileName<<endl;
    ifstream fin;
    fin.open(aFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aFileName);
    ProgressBar progress_bar;
    while (!fin.eof() && fin.good()){
      string line;
      getline(fin,line);
      stringstream ss;
      ss<<line<<endl;
      while (!ss.eof() && ss.good()){
	string target;
	ss>>target;
	if (target!=""){	    
	  oStringList.push_back(target);
	  progress_bar.Count();
	}
      }
    }
    fin.close();
  }

  void InputAndOutput(const string& aInputFileName){
    Load(aInputFileName,true);
  }

  void Input(const string& aInputFileName){
    Load(aInputFileName,false);
  }

  void Load(const string& aInputFileName, bool aDirectProcess){
    string ofname=aInputFileName+".feature";
    ofstream ofs_f(ofname.c_str());
    ofname=aInputFileName+".feature_bin";
    ofstream ofs_fb(ofname.c_str());

    cout<<"Reading gspan data and computing features"<<endl;
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    ProgressBar progress_bar;
    while (!fin.eof()){
      GraphClass G;
      SetGraphFromFileGSPAN(fin, G, mMode);
      SVector x;
      Generate(G,x);
      if (aDirectProcess){
	ofs_f<<x;
	x.save(ofs_fb);
      } else { 
	mDataset.push_back(x);
      }
      progress_bar.Count();
    }
    fin.close();
  }

  void InputSparse(const string& aInputFileName, string aMode){
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    InputSparse(fin,aMode,mDataset);
    fin.close();
  }

  void InputSparse(const string& aInputFileName, string aMode, vector<SVector>& oDataset){
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    InputSparse(fin,aMode,oDataset);
    fin.close();
  }

  void InputSparse(ifstream& aFin, string aMode, vector<SVector>& oDataset){
    cout<<"Reading file in "<<aMode<<" mode"<<endl;
    ProgressBar progress_bar;
    while (!aFin.eof() && aFin.good()){
      SVector x;
      if (aMode=="binary") x.load(aFin);
      else ParseASCIILine2Vector(aFin,x);
      if (InstanceIsValid(x)==true) {
	oDataset.push_back(x);
	progress_bar.Count();
      } else {}//discard non valid instances
    }
  }

  inline void ParseASCIILine2Vector(ifstream& aFin, SVector& aX){
    string line;
    getline(aFin,line);
    if (line=="") return;
    stringstream ss;
    ss<<line<<endl;
    while (!ss.eof() && ss.good()){
      string key_value;
      ss>>key_value;
      size_t limit=key_value.find_first_of(":", 0);
      string key=key_value.substr(0, limit);
      string value=key_value.substr(limit+1, key_value.size());
      unsigned key_int=stream_cast<unsigned>(key);
      double val_real=stream_cast<double>(value);
      aX.set(key_int,val_real);
    }
  }

  inline  bool InstanceIsValid(SVector& aX){
    bool is_valid=false;
    //if (dot(aX,aX)>0) is_valid=true;
    if (aX.sparse_size()>0) is_valid=true;
    return is_valid;
  }

  void SetGraphFromFileGSPAN(istream& in, GraphClass& oG, string aMode=""){
    map<string,int> index_map_nominal_to_real;
    string line;
    getline(in,line);
    assert(line[0]=='t');//first line must have as first char a 't'
    while(!in.eof() && in.good() && in.peek()!='t' && getline(in,line)){//read until next 't' or end of file
      stringstream ss;
      ss<<line<<endl;
      char code;
      ss>>code;
      if (code=='v'){
	//extract vertex id and make map nominal_id -> real_id
	string nominal_vertex_index;
	ss>>nominal_vertex_index;
	unsigned real_vertex_index=oG.InsertVertex();
	index_map_nominal_to_real[nominal_vertex_index]=real_vertex_index;
	//label
	vector<string> vertex_symbolic_attribute_list;
	string label;
	ss>>label;
	////vertex_symbolic_attribute_list.push_back("nt");
	vertex_symbolic_attribute_list.push_back(label);
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	if (aMode=="3D"){
	  double x,y,z;
	  ss>>x>>y>>z;
	  vector<double> vertex_numeric_attribute_list;
	  vertex_numeric_attribute_list.push_back(x);
	  vertex_numeric_attribute_list.push_back(y);
	  vertex_numeric_attribute_list.push_back(z);
	  oG.SetVertexNumericAttributeList(real_vertex_index, vertex_numeric_attribute_list);
	}

	//status
	vector<bool> status; 
	 
	status.push_back(true);//kernel point
	status.push_back(true);//kind
	if (aMode=="viewpoint"){//if mode=viewpoint then assume extra vertex label is \in {0,1} representing viewpoint flag
	  bool vp;
	  ss>>vp;
	  status.push_back(vp);//viewpoint
	}	  
	else status.push_back(true);//viewpoint
	status.push_back(false);//dead
	status.push_back(false);//forbidden distance
	status.push_back(false);//forbidden neighborhood
	oG.SetVertexStatusAttributeList(real_vertex_index, status);
      } else if (code=='e'){
	//extract src and dest vertex id 
	string nominal_src_index,nominal_dest_index;
	string label;
	double weight;
	ss>>nominal_src_index>>nominal_dest_index>>label>>weight;
	unsigned real_src_index=index_map_nominal_to_real[nominal_src_index];
	unsigned real_dest_index=index_map_nominal_to_real[nominal_dest_index];
	unsigned edge_index=oG.InsertEdge(real_src_index,real_dest_index);
	unsigned reverse_edge_index=oG.InsertEdge(real_dest_index,real_src_index);
	vector<string> edge_symbolic_attribute_list;
	edge_symbolic_attribute_list.push_back(label);
	vector<double> edge_numeric_attribute_list;
	edge_numeric_attribute_list.push_back(weight);
	oG.SetEdgeSymbolicAttributeList(edge_index,edge_symbolic_attribute_list);
	oG.SetEdgeSymbolicAttributeList(reverse_edge_index,edge_symbolic_attribute_list);
	oG.SetEdgeNumericAttributeList(edge_index,edge_numeric_attribute_list);
	oG.SetEdgeNumericAttributeList(reverse_edge_index,edge_numeric_attribute_list);
      } else {}//NOTE: ignore other markers
    }
  }

  void ComputeBinDataStructure(unsigned aNumHashFunctions){
    string ofname="hash_encoding";
    ofstream of(ofname.c_str());
        
    cout<<"Computing bin data structure..."<<endl;
    ProgressBar progress_bar;

    //init structure
    mBinDataStructure.clear();
    for (unsigned k=0;k<aNumHashFunctions;++k)
      mBinDataStructure.push_back(multimap<unsigned,unsigned>());
    
    //fill structure
    for (unsigned i=0;i<mDataset.size();++i){
      vector<unsigned> min_list=ComputeHashSignature(i,aNumHashFunctions);
      if (mOutputHashEncoding) {for(unsigned j=0;j<min_list.size();j++) of<<min_list[j]<<" ";of<<endl;}
      for (unsigned k=0;k<aNumHashFunctions;++k)
	mBinDataStructure[k].insert(make_pair(min_list[k],i));
      progress_bar.Count();
    }
  }

  inline  vector<unsigned> ComputeHashSignature(unsigned aID, unsigned aNumHashFunctions){
    if (mSignatureMap.count(aID)>0) return mSignatureMap[aID];
    else { 
      vector<unsigned> signature=ComputeHashSignature(mDataset[aID],aNumHashFunctions);
      mSignatureMap[aID]=signature;
      return signature;
    }
  }

  inline  vector<unsigned> ComputeHashSignature(SVector& aX, unsigned aNumHashFunctions){
    unsigned effective_num_hash_functions=aNumHashFunctions*mHashFactor;
    const double A=sqrt(2)-1;
    const unsigned MAXUNSIGNED=2<<30;
    vector<unsigned> signature;
    //special case for first hash function just consider first feature id
    //signature.push_back(mDataset[aID].extract_component(0).first);
    //for all other hash functions rehash and select the minimum
    for (unsigned k=0;k<effective_num_hash_functions;++k)
      signature.push_back(MAXUNSIGNED);
    unsigned size=(unsigned)aX.sparse_size();
    for (unsigned f=0;f<size;++f){    
      unsigned hash_id=aX.extract_component(f).first;
      for (unsigned k=0;k<effective_num_hash_functions;++k){
	unsigned new_hash=IntHash(hash_id*(k+1)*A);
	if (signature[k]>new_hash) signature[k]=new_hash;
      }
    }
    //compact signature
    vector<unsigned> compact_signature;
    for (unsigned i=0;i<signature.size();i=i+mHashFactor){
      unsigned new_hash=0;
      for(unsigned j=0;j<mHashFactor;j++)
	new_hash+=signature[i+j];
      compact_signature.push_back(new_hash);
    }
    return compact_signature;
  }

  void ComputeInvertedFeatureInstanceIndex(){
    cout<<"Computing inverted data structure..."<<endl;
    ProgressBar progress_bar;
    for (unsigned i=0;i<mDataset.size();++i){
      for (unsigned k=0;k<(unsigned)mDataset[i].sparse_size();++k){
	unsigned hash_id=mDataset[i].extract_component(k).first;
	mInvertedIndex.insert(make_pair(hash_id,i));
      }
      progress_bar.Count();
    }
  }

  vector<unsigned> ComputeFastNeighbourhood(unsigned aID){
    map<unsigned,int> neighborhood;
    for (unsigned k=0;k<(unsigned)mDataset[aID].sparse_size();++k){
      unsigned hash_id=mDataset[aID].extract_component(k).first;
      //return equal range wrt hash_id
      pair<multimap<unsigned,unsigned>::iterator,multimap<unsigned,unsigned>::iterator> erange=mInvertedIndex.equal_range(hash_id);
      //insert into neighborhood
      for (multimap<unsigned,unsigned>::iterator it=erange.first;it!=erange.second;++it){
	unsigned instance_id=it->second;
	//count number of occurences
	if (neighborhood.count(instance_id)>0) neighborhood[instance_id]++;
	else neighborhood[instance_id]=1;
      }
    }
    return TrimNeighbourhood(neighborhood);
  }

  vector<unsigned> ComputeApproximateNeighbourhood(unsigned aID){
    vector<unsigned> hash_signature=ComputeHashSignature(aID,mNumHashFunctions);
    return ComputeApproximateNeighbourhood(hash_signature);
  }

  vector<unsigned> ComputeApproximateNeighbourhood(const vector<unsigned>& aInstance){
    map<unsigned,int> neighborhood;
    vector<pair<unsigned, double> > vec;
    for (unsigned k=0;k<mNumHashFunctions;++k){
      unsigned hash_id=aInstance[k];
      //return equal range wrt hash_id
      pair<multimap<unsigned,unsigned>::iterator,multimap<unsigned,unsigned>::iterator> erange=mBinDataStructure[k].equal_range(hash_id);
      //fill neighborhood set counting number of occurrences
      for (multimap<unsigned,unsigned>::iterator it=erange.first;it!=erange.second;++it){
	unsigned instance_id=it->second;
	if (neighborhood.count(instance_id)>0) neighborhood[instance_id]++;
	else neighborhood[instance_id]=1;
      }
    }
    return TrimNeighbourhood(neighborhood);
  }

  vector<unsigned> ComputeTrueNeighbourhood(const SVector& aX, unsigned aSize){
    vector<pair<double, unsigned> > rank_list;
    for (unsigned i=0;i<mDataset.size();++i){
      double k=dot(aX,mDataset[i]);
      rank_list.push_back(make_pair(-k,i));
    }
    unsigned effective_size=min((unsigned)rank_list.size(),aSize);
    partial_sort(rank_list.begin(),rank_list.begin()+effective_size,rank_list.end());
    vector<unsigned> neighbor;
    for (unsigned j=0;j<effective_size;j++)
      neighbor.push_back(rank_list[j].second);
    return neighbor;
  }

  vector<unsigned> TrimNeighbourhood(map<unsigned,int>& aNeighborhood){
    //given a list of neighbours with an associated occurences count, return only a fraction of the highest count ones
    vector<unsigned> neighborhood_list;
    if (mEccessNeighbourSizeFactor>0){
      //sort by num occurences
      vector<pair<int,unsigned> > count_list;
      for (map<unsigned,int>::const_iterator it=aNeighborhood.begin();it!=aNeighborhood.end();++it){
	unsigned id=it->first;
	int count=it->second;
	count_list.push_back(make_pair(-count,id));//NOTE:-count to sort from highest to lowest
      }
      unsigned effective_size=min((unsigned)count_list.size(),(unsigned)(mEccessNeighbourSizeFactor*mNeighborhoodSize));
      //partial_sort(count_list.begin(),count_list.begin()+effective_size,count_list.end());
      sort(count_list.begin(),count_list.end());
      for (unsigned i=0;i<effective_size;++i)
	neighborhood_list.push_back(count_list[i].second);
    } else {//if mEccessNeighbourSizeFactor==0 then just return all the ids in the  neighborhood 
      for (map<unsigned,int>::const_iterator it=aNeighborhood.begin();it!=aNeighborhood.end();++it){
	neighborhood_list.push_back(it->first);
      }
    }
    return neighborhood_list;
  }

  double ComputeApproximateDensity(unsigned aID, unsigned aSize){
    vector<unsigned> approximate_neighborhood=ComputeApproximateNeighbourhood(aID);
    vector<pair<double,unsigned> > sim_list;
    //compute kernel pairs between aID and all elements in aApproximateNeighborhood
    for (unsigned j=0;j<approximate_neighborhood.size();j++){
      unsigned u=aID;
      unsigned v=approximate_neighborhood[j];
      if (u!=v){
	double k_uv=Kernel(u,v);
	sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
      }
    }
    //sort and take truly most similar 
    sort(sim_list.begin(),sim_list.end());
    if ((unsigned)sim_list.size()<aSize) return 0;//NOTE: use only instances for which there are approximate neighborhoods of size=aSize 
    unsigned effective_size=min(aSize,(unsigned)sim_list.size());
    double density=0;
    for (unsigned k=0;k<effective_size;++k)
      density+=-sim_list[k].first;//note: restore the similarity as a positive numberx
    double relative_density=effective_size==0?0:density/effective_size;
    //cout<<aID<<" Size: "<<aSize<<" ApproxNeighborhoodSize: "<<approximate_neighborhood.size()<<" EstimatedDensity: "<<relative_density<<endl;////////////////////////////
    return relative_density;
  }

  double ComputeTrueDensity(unsigned aID){
    vector<pair<double,unsigned> > sim_list;
    double density=0;
    for (unsigned j=0;j<mDataset.size();j++){
      if (aID!=j){
	double k_ij=Kernel(aID,j);
	density+=k_ij;
      }
    }
    return density/mDataset.size();
  }

  vector<unsigned> ComputeRandomSampleHighDensityInstances(unsigned aNeighborhoodSize, unsigned aSampleSize){
    //compute density for each instance in dataset
    cout<<"Computing approximate density information for each instance."<<endl;////
    vector<double> density_list;
    double max_density_val=ComputeApproximateDensity(0,aNeighborhoodSize);
    unsigned max_density_id=0;
    double min_density_val=max_density_val;
    unsigned min_density_id=0;
    unsigned non_zero_density_count=0;
    {
      ProgressBar progress_bar;
      progress_bar.Count();
      for (unsigned i=1;i<mDataset.size();++i){
	double density=ComputeApproximateDensity(i,aNeighborhoodSize);
	density_list.push_back(density);
	progress_bar.Count();
	if (density>0){
	  non_zero_density_count++;
	  if (max_density_val<density) {
	    max_density_val=density;
	    max_density_id=i;
	  }
	  if (min_density_val>density) {
	    min_density_val=density;
	    min_density_id=i;
	  }
	}
      }
    }

    //sample tentative centers
    cout<<"Sampling high density instances."<<endl;////
    set<unsigned> center_sample_set;
    center_sample_set.insert(max_density_id);
    const unsigned SAMPLE_SIZE_FACTOR=4;//extract a set of tentative high density instances of size SAMPLE_SIZE_FACTOR*aSampleSize 
    const int MAX_NUM_DATASET_SCANS=5000;
    int counter=MAX_NUM_DATASET_SCANS;
    unsigned effective_size_center_sample=min(non_zero_density_count,aSampleSize*SAMPLE_SIZE_FACTOR);
    {
      ProgressBar progress_bar;
      while (counter>0 && center_sample_set.size()<effective_size_center_sample){
	for (unsigned i=0;i<density_list.size();++i){
	  double d=(density_list[i]-min_density_val)/(max_density_val-min_density_val);//Normalize density between 1=max density and 0=min density
	  double r=random01();
	  if (d > r) {//bias current instance selection so that it is not too similar to previously selected instances
	    double avg_sim=ComputeAverageSimilarityToSet(i,center_sample_set);
	    double sim=pow(avg_sim,mNonRedundantFilter);
	    double r=random01();
	    if ( sim < r ){
	      center_sample_set.insert(i);
	      density_list[i]=0;
	      progress_bar.Count();
	    } 
	  }
	}
	counter--;
      }
    }

    //compute the true density for the selected sample and retain only the highest density instances
    cout<<"Compute true density for selected "<<center_sample_set.size()<<" centers."<<endl;////
    vector<pair<double,unsigned> > sim_list;
    {
      ProgressBar progress_bar;
      for (set<unsigned>::iterator it=center_sample_set.begin();it!=center_sample_set.end();++it){
	unsigned id=(*it);
	double value=-ComputeTrueDensity(id);//NOTE:-density to sort starting from highest true density
	sim_list.push_back(make_pair(value,id));
	progress_bar.Count();
      }
    }
    sort(sim_list.begin(),sim_list.end());
    //sample one every SAMPLE_SIZE_FACTOR entries to guarantee centers diversity
    vector<unsigned> sample;
    for (unsigned k=0;k<sim_list.size() && sample.size()<aSampleSize;++k)
      if (k%SAMPLE_SIZE_FACTOR==0) sample.push_back(sim_list[k].second);
    return sample;
  }

  double ComputeAverageSimilarityToSet(unsigned aID, set<unsigned>& aSet){
    double avg_sim=0;
    for (set<unsigned>::iterator it=aSet.begin();it!=aSet.end();++it){
      unsigned id=(*it);
      avg_sim+=Kernel(aID,id);
    }
    return avg_sim/aSet.size();
  }

  vector<unsigned> ComputeNeighborhoodRanking(unsigned aID){
    vector<pair<double,unsigned> > sim_list;
    for (unsigned i=0;i<mDataset.size();++i){
      if (i!=aID) {
	double k=Kernel(aID,i);
	sim_list.push_back(make_pair(-k,i));//note: use -k to sort in decreasing order
      }
    }
    //sort and take most similar <effective_size>
    sort(sim_list.begin(),sim_list.end());
    vector<unsigned> neighborhood;
    for (unsigned t=0;t<sim_list.size();++t)
      neighborhood.push_back(sim_list[t].second);
    return neighborhood;
  }

  void OutputCluster(ostream& out){
    vector<unsigned> density_center_list=ComputeRandomSampleHighDensityInstances(mNeighborhoodSize,mSampleSize);
    cout<<"Compute neighbourhood for selected "<<density_center_list.size()<<" cluster centers."<<endl;////
    ProgressBar progress_bar;
    for (unsigned i=0;i<density_center_list.size();++i){
      unsigned id=density_center_list[i];
      vector<unsigned> neighborhood=ComputeNeighborhoodRanking(id);
      progress_bar.Count();
      out<<id+1<<"   ";      //note: output numbering starts from 1 (not from 0)
      for (unsigned i=0;i<neighborhood.size(); ++i )
	out<<neighborhood[i]+1<<" ";       //note: output numbering starts from 1 (not from 0)
      out<<endl;
    }
  }

  void OutputApproximateKNN(ostream& out){
    cout<<"Compute approximate nearest neighbours."<<endl;////
    ProgressBar progress_bar;
    for (unsigned u=0;u<mDataset.size();++u){
      vector<unsigned> approximate_neighborhood=ComputeApproximateNeighbourhood(u);

      if (mTrueSort==true){//sort all instances in neighbourhood according to true similarity and take only the mNeighborhoodSize nearest ones
	vector<pair<double,unsigned> > sim_list;
	//compute kernel pairs between aID and all elements in aApproximateNeighborhood
	for (unsigned j=0;j<approximate_neighborhood.size();j++){
	  unsigned v=approximate_neighborhood[j];
	  double k_uv=Kernel(u,v);
	  sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
	}
	//sort and take truly most similar 
	sort(sim_list.begin(),sim_list.end());
	unsigned effective_size=min(mNeighborhoodSize,(unsigned)sim_list.size());
	for (unsigned k=0;k<effective_size;++k){
	  out<<sim_list[k].second+1<<" ";//NOTE: numbering starts from 1
	}
	out<<endl;
      } else {//do not sort neighbours according to true similarity, just output them
	for (unsigned t=0;t<approximate_neighborhood.size();++t){
	  out<<approximate_neighborhood[t]+1<<" ";//NOTE: numbering starts from 1
	}
	out<<endl;
      }
      progress_bar.Count();
    }
  }

  void OutputFastKNN(ostream& out){
    cout<<"Compute fast nearest neighbours."<<endl;////
    ProgressBar progress_bar;
    for (unsigned u=0;u<mDataset.size();++u){
      vector<unsigned> neighborhood=ComputeFastNeighbourhood(u);
      vector<pair<double,unsigned> > sim_list;
      //compute kernel pairs between aID and all elements in neighborhood
      for (unsigned j=0;j<neighborhood.size();j++){
	unsigned v=neighborhood[j];
	double k_uv=Kernel(u,v);
	sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
      }
      //sort and take truly most similar 
      sort(sim_list.begin(),sim_list.end());
      unsigned effective_size=min(mNeighborhoodSize,(unsigned)sim_list.size());
      for (unsigned k=0;k<effective_size;++k){
	out<<sim_list[k].second+1<<" ";//NOTE: numbering starts from 1
      }
      out<<endl;
      progress_bar.Count();
    }
  }

  void OutputTrueKNN(ostream& out){
    cout<<"Compute true nearest neighbours."<<endl;////
    ProgressBar progress_bar;
    for (unsigned u=0;u<mDataset.size();++u){
      vector<pair<double,unsigned> > sim_list;
      //compute kernel pairs between aID and all elements 
      for (unsigned v=0;v<mDataset.size();v++){
	double k_uv=Kernel(u,v);
	sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
      }
      //sort and take truly most similar 
      sort(sim_list.begin(),sim_list.end());
      for (unsigned k=0;k<mNeighborhoodSize;++k){
	out<<sim_list[k].second+1<<" ";//NOTE: numbering starts from 1
      }
      out<<endl;
      progress_bar.Count();
    }
  }

  void OutputApproximateKNNPrediction(ostream& out, string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName){
    cout<<"Compute approximate nearest neighbour prediction of test instances."<<endl;////
    ProgressBar progress_bar;

    //read test instances
    vector<SVector> test_dataset;
    InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

    //read train targets
    vector<string> train_target_list;
    InputStringList(aTrainTargetInputFileName,train_target_list);

    cout<<"Computing k-NN predictions."<<endl;
    //for each test instance
    for (unsigned u=0;u<test_dataset.size();++u){
      //extract signature
      vector<unsigned> hash_signature=ComputeHashSignature(test_dataset[u],mNumHashFunctions);
      //extract knn
      vector<unsigned> approximate_neighborhood=ComputeApproximateNeighbourhood(hash_signature);
      string prediction=KNNPredict(approximate_neighborhood,train_target_list);
      out<<prediction<<endl;
      progress_bar.Count();
    }
  }

  void OutputTrueKNNPrediction(ostream& out, string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName){
    cout<<"Compute true nearest neighbour prediction of test instances."<<endl;////
    ProgressBar progress_bar;

    //read test instances
    vector<SVector> test_dataset;
    InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

    //read train targets
    vector<string> train_target_list;
    InputStringList(aTrainTargetInputFileName,train_target_list);

    cout<<"Computing k-NN predictions on "<<test_dataset.size()<<" instances."<<endl;
    //for each test instance
    for (unsigned u=0;u<test_dataset.size();++u){
      //extract knn
      vector<unsigned> neighborhood=ComputeTrueNeighbourhood(test_dataset[u],mNeighborhoodSize);
      //compute majority class
      string prediction=KNNPredict(neighborhood,train_target_list);
      out<<prediction<<endl;
      progress_bar.Count();
    }
  }

  string KNNPredict(const vector<unsigned>& aNeighborhood, const vector<string>& aTargetList)const{
    //compute histogram of targets in neighborhood
    map<string,unsigned> histogram;
    for (unsigned i=0;i<aNeighborhood.size();++i){
      unsigned nn_id=aNeighborhood[i];
      assert(nn_id<aTargetList.size());      
      string predicted_target=aTargetList[nn_id];
      if (histogram.count(predicted_target)==0) histogram[predicted_target]=1;
      else histogram[predicted_target]++;
    }
    //compute majority vote for target
    string max_target=aTargetList[0];//initialization with one arbitrary target
    unsigned max_val=0;
    for (map<string,unsigned>::const_iterator it=histogram.begin();it!=histogram.end();++it){
      string target=it->first;
      unsigned vote=it->second;
      if (max_val<vote) {
	max_val=vote;
	max_target=target;
      }
    }
    return max_target;
  }

  void OutputDatasetSplit(string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName, unsigned aNumCenters, unsigned aSizeThreshold, double aImbalanceTolerance){
    cout<<"Compute train split."<<endl;////
    //read test instances
    vector<SVector> test_dataset;
    InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

    //read train targets
    vector<string> train_target_list;
    InputStringList(aTrainTargetInputFileName,train_target_list);
    
    //randomly shuffle test instances
    vector<unsigned> test_id_list(test_dataset.size());
    for (unsigned i=0;i<test_dataset.size(); ++i ) test_id_list.push_back(i);
    random_shuffle(test_id_list.begin(),test_id_list.end());
    //iterate until the prespecified number of centers that satisfy the required considtions are found
    vector<vector<double> > valid_center_list;
    cout<<"Computing test centers."<<endl;
    ProgressBar progress_bar;
    for (unsigned i=0; i<test_id_list.size();++i){
      progress_bar.Count();
      //check if it is possible to extract a valid train subset from the selected test instance
      unsigned test_id=test_id_list[i];
      //sort all training set according to distance from selected test intance
      vector<unsigned> neighborhood=ComputeTrueNeighbourhood(test_dataset[test_id],train_target_list.size());
      //compute the sorted sequence of targets
      vector<string> neighborhood_target;
      for (unsigned j=0;j<neighborhood.size();j++)
	neighborhood_target.push_back(train_target_list[neighborhood[j]]);
      //return the largest number of neighbour such that the target distribution is uniform
      double positive_class_count=0;
      double negative_class_count=0;
      unsigned ksize=0;
      double kvalue=0;
      double nbalance=0;
      for (unsigned j=0;j<neighborhood_target.size();j++){
	if (neighborhood_target[j]=="1") positive_class_count++;
	else if (neighborhood_target[j]=="-1") negative_class_count++;
	double balance=(negative_class_count>0)?positive_class_count/negative_class_count:0;
	if (balance>1/aImbalanceTolerance && balance<aImbalanceTolerance) {
	  ksize=j;
	  unsigned train_id=neighborhood[j];
	  kvalue=dot(mDataset[train_id],test_dataset[test_id]);
	  nbalance=balance;
	}
      }
      //if a valid neighbourhood has been found then store it
      if (ksize>aSizeThreshold) {
	vector<double> item;
	item.push_back((double)test_id);
	item.push_back((double)ksize);
	item.push_back(kvalue);
	item.push_back(nbalance);
	valid_center_list.push_back(item);
      }
    }

    //iterate over all aNumCenters valid centers
    unsigned dataset_counter=0;
    for (unsigned i=0; i<valid_center_list.size() && dataset_counter<aNumCenters;++i){
      assert(valid_center_list[i].size()==4);
      unsigned center_id=(unsigned)valid_center_list[i][0];
      //unsigned ksize=(unsigned)valid_center_list[i][1];
      double kvalue=valid_center_list[i][2];
      double nbalance=valid_center_list[i][3];

      //iterate over test set and collect test instances that are within kvalue from the center 
      vector<unsigned> part_test_list;
      for (unsigned j=0;j<test_dataset.size();++j){
	double k=dot(test_dataset[j],test_dataset[center_id]);
	if (k>=kvalue) part_test_list.push_back(j);
      }
      //if condistions are met, then output train-test datsets
      if (part_test_list.size()>aSizeThreshold){
	//output test set
	string ofnamets=aSparseASCIITestInputFileName+".test_"+stream_cast<string>(dataset_counter);
	ofstream ofsts(ofnamets.c_str());
	for (unsigned t=0;t<part_test_list.size();++t) ofsts<<part_test_list[t]+1<<endl;//+1 as ids range from 1 and not from 0
	ofsts.close();

	//collect training instances that are in the balanced neighborhood
	vector<unsigned> part_train_list;
	for (unsigned j=0;j<mDataset.size();++j){
	  double k=dot(mDataset[j],test_dataset[center_id]);
	  if (k>=kvalue) part_train_list.push_back(j);
	}
	//output training set
	string ofnametr=aSparseASCIITestInputFileName+".train_"+stream_cast<string>(dataset_counter);
	ofstream ofstr(ofnametr.c_str());
	for (unsigned t=0;t<part_train_list.size();++t) {
	  
	  ofstr<<part_train_list[t]+1<<endl;//+1 as ids range from 1 and not from 0
	}
	ofstr.close();

	//logging info
	LOGF<<dataset_counter<<" center_id:"<<center_id+1<<" kval:"<<kvalue<<" #train:"<<part_train_list.size()<<" #test:"<<part_test_list.size()<<" balance:"<<nbalance<<endl;
      
	dataset_counter++;
      }
    }
    cout<<endl<<"Computed "<<dataset_counter<<" valid centers"<<endl;
  }

  
  void OutputKernel(ostream& out){
    cout<<"Compute kernel matrix."<<endl;////
    ProgressBar progress_bar;
    for (unsigned i=0;i<mDataset.size(); ++i ){
      for (unsigned j=0;j<mDataset.size();j++)
	out<<Kernel(i,j)<<" ";
      out<<endl;
      progress_bar.Count();
    }
  }

  void Output(ostream& out){
    for (unsigned i=0;i<mDataset.size(); ++i )
      out<<mDataset[i];
  }

  void OutputFeatureMap(string aFileName)const{
    {
      string ofname=aFileName+".feature_map";
      ofstream of(ofname.c_str());
      pmFeatureGenerator->OutputFeatureMap(of);
    }
    {
      string ofname=aFileName+".soft_feature_map";
      ofstream of(ofname.c_str());
      pmSoftFeatureGenerator->OutputFeatureMap(of);
    }
  }

  double Kernel(unsigned aI, unsigned aJ){
    if (mCaching){
      unsigned i=min(aI,aJ);
      unsigned j=max(aI,aJ);
      pair<unsigned,unsigned> key=make_pair(i,j); 
      if (mKernelMap.count(key)==0){
	double value=dot(mDataset[i],mDataset[j]);
	mKernelMap[key]=value;
      }
      return mKernelMap[key];
    } else 
      return dot(mDataset[aI],mDataset[aJ]);
  }
};

//---------------------------------------------------------------------------------
int main(int argc, char** argv){
  Timer T;
  srand(time(NULL));
  LOGF<<"--------------------------------------------------------------------------------"<<endl;
  LOGF<<CREDITS<<endl;
  time_t rawtime_start; time ( &rawtime_start ); LOGF<<"Start logging: "<<asctime ( localtime ( &rawtime_start ) )<<endl;
  LOGF<<"Command line: ";for(int i=0;i<argc; ++i ) LOGF<<stream_cast<string>(argv[i])<<" ";LOGF<<endl;
  try
    {
      ParameterWrapperClass P;
      P.Init(argc,argv);

      string mode="";
      //factory
      NSPDK_FeatureGenerator fg("nspdk");
      ANSPDK_FeatureGenerator afg("anspdk");
      ALNSPDK_FeatureGenerator alfg("alnspdk");
      MNSPDK_FeatureGenerator mfg("mnspdk");
      NSPDK3D_FeatureGenerator fg3d("nspdk3d");

      NSPDK_FeatureGenerator* pfg;
      if (P.mType=="nspdk") pfg=&fg;
      else if (P.mType=="anspdk") pfg=&afg;
      else if (P.mType=="alnspdk") pfg=&alfg;
      else if (P.mType=="mnspdk") pfg=&mfg;
      else if (P.mType=="nspdk3d") {
	pfg=&fg3d;
	mode="3D";
      }
      else if (P.mType=="nspdkvp") {
	pfg=&fg;
	mode="viewpoint";
      }
      else throw range_error("Unknown feature generator type:"+P.mType);

      pfg->set_flag("radius",stream_cast<string>(P.mRadiusMax));
      pfg->set_flag("distance",stream_cast<string>(P.mDistanceMax));
      pfg->set_flag("match_type",stream_cast<string>(P.mMatchingType));
      pfg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
      pfg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
      pfg->set_flag("verbosity",stream_cast<string>(P.mDebug));
      if (P.mMinKernel) pfg->set_flag("min_kernel","true");
      if (!P.mNormalization) pfg->set_flag("normalization","false");

      //feature generator for additional soft mode
      NSPDK_FeatureGenerator soft_fg("soft_nspdk");
      NSPDK_FeatureGenerator* p_soft_fg=&soft_fg;
      p_soft_fg->set_flag("radius",stream_cast<string>(P.mRadiusMax));
      p_soft_fg->set_flag("distance",stream_cast<string>(P.mDistanceMax));
      p_soft_fg->set_flag("match_type","soft");
      p_soft_fg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
      p_soft_fg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
      p_soft_fg->set_flag("verbosity",stream_cast<string>(P.mDebug));
      if (P.mMinKernel) p_soft_fg->set_flag("min_kernel","true");
      if (!P.mNormalization) p_soft_fg->set_flag("normalization","false");


      string ofname;
      //main process
      MainProcessClass C(pfg,
			 p_soft_fg,
			 P.mNumHashFunctions,
			 P.mHashFactor,
			 P.mNumNearestNeighbors,
			 P.mSampleSize,
			 P.mNonRedundantFilter,
			 P.mAbstractFeatures,
			 mode, 
			 P.mCaching, 
			 P.mTrueSort, 
			 P.mEccessNeighbourSizeFactor,
			 P.mMixedMatchingType,
			 P.mOutputHashEncoding);
      if (P.mOutputFeatures) {
	C.InputAndOutput(P.mGspanInputFileName);
	if (P.mOutputFeatureMap)
	  C.OutputFeatureMap(P.mGspanInputFileName);
      } else {
	if (P.mSparseBinaryInputFileName!="")
	  C.InputSparse(P.mSparseBinaryInputFileName,"binary");
	else if (P.mSparseASCIIInputFileName!="")
	  C.InputSparse(P.mSparseASCIIInputFileName,"ascii");
	else if (P.mGspanInputFileName!="")
	  C.Input(P.mGspanInputFileName);
	else
	  throw range_error("ERROR:No input file name specified");

	if (P.mOutputCluster){
	  C.ComputeBinDataStructure(P.mNumHashFunctions);
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".fast_cluster";
	  ofstream ofs_fc(ofname.c_str());
	  C.OutputCluster(ofs_fc);
	}

	if (P.mOutputApproximateKNN){
	  C.ComputeBinDataStructure(P.mNumHashFunctions);
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".approx_knn";
	  ofstream ofs_aknn(ofname.c_str());
	  C.OutputApproximateKNN(ofs_aknn);
	}

	if (P.mOutputFastKNN){
	  C.ComputeInvertedFeatureInstanceIndex();
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".fast_knn";
	  ofstream ofs_fknn(ofname.c_str());
	  C.OutputFastKNN(ofs_fknn);
	}

	if (P.mOutputTrueKNN){
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".knn";
	  ofstream ofs_fknn(ofname.c_str());
	  C.OutputTrueKNN(ofs_fknn);
	}

	if (P.mOutputKernel){
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".kernel";
	  ofstream ofs_fk(ofname.c_str());
	  C.OutputKernel(ofs_fk);
	}

	if (P.mOutputApproximateKNNPrediction){
	  C.ComputeBinDataStructure(P.mNumHashFunctions);
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".approx_knn_prediction";
	  ofstream ofs_knnp(ofname.c_str());
	  C.OutputApproximateKNNPrediction(ofs_knnp,P.mSparseASCIITestInputFileName,P.mTrainTargetInputFileName);
	}

	if (P.mOutputTrueKNNPrediction){
	  ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".knn_prediction";
	  ofstream ofs_knnp(ofname.c_str());
	  C.OutputTrueKNNPrediction(ofs_knnp,P.mSparseASCIITestInputFileName,P.mTrainTargetInputFileName);
	}
	if (P.mOutputDatasetSplit){
	  C.OutputDatasetSplit(P.mSparseASCIITestInputFileName,P.mTrainTargetInputFileName,P.mNumCenters, P.mSizeThreshold,P.mImbalanceTolerance);
	}
      }
    }
  catch(exception& e)
    {
      cerr<<e.what();
      LOGF<<e.what()<<endl;
    }
  time_t rawtime_end; time ( &rawtime_end );LOGF<<"End logging: "<<asctime ( localtime ( &rawtime_end ) )<<endl;
  LOGF<<"Time elapsed: "<<T.end()<<endl;
  return 0;
}
