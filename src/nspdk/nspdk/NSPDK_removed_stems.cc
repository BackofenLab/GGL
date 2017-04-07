#include "nspdk/Utility.h"
#include "nspdk/SVector.h"
#include "nspdk/Histogram.h"
#include "nspdk/BaseGraphClass.h"
#include "nspdk/GraphClass.h"
#include "nspdk/NSPDK_FeatureGenerator.h"



namespace nspdk {

using namespace std;


//---------------------------------------------------------------------------------
const string CREDITS("Name: NeighborhoodSubgraphPaiwiseDistanceKernel\nVersion: 6.1\nProgrammer: Fabrizio Costa\nLast update: July 21 2010");
//---------------------------------------------------------------------------------
ofstream LOGF("log",std::ios_base::app);

//FlagsService& The_FlagsService= FlagsService::get_instance();

class ParameterWrapperClass{
public:
  ParameterWrapperClass():mInputFileName(""),
			  mRadiusMax(1),
			  mRadiusMin(0),
			  mDistanceMax(4),
			  mDistanceMin(0),
			  mMatchingType("hard"),
			  mVerbose(false),
			  mFeatureBitSize(30),
			  mMinKernel(false),
			  mType("nspdk"),
			  mOutputKernel(false),  
			  mNormalization(true),
			  mNumNearestNeighbors(7),
			  mSaveDot(false),
			  mAlpha(.5),
			  mZThresholdLink(1),
			  mZThresholdCentrality(1),
			  mZThresholdSubtreeSize(1),
			  mSizeNeighborhood(7),
			  mRootDistanceLimit(3),
			  mMarkParts(false),
			  mOutputFeatures(false),
			  mDebug(false){}
  
  void Usage(string aCommandName){
    cerr<<"Usage: "<<aCommandName<<endl
	<<"[-f <file name gspan format> for input from file (default: stdin)]"<<endl
	<<"[-of flag to output only feature encoding (default: false -> output similarity matrix)]"<<endl
	<<"[-knn <num nearest neighbors> (default: 7)]"<<endl
	<<"[-a <alpha value to tradeoff centrality with dissimilarity> (default: .5)]"<<endl
	<<"[-zl <Z threshold for Link> (default: 1)]"<<endl
	<<"[-zc <Z threshold for Centrality> (default: 1)]"<<endl
	<<"[-zs <Z threshold for subtree size> (default: 1)]"<<endl
	<<"[-sn <size neighborhood to estimate average pairwise similarity reference (default: 7)]"<<endl
	<<"[-rdl <root distance limit (default: 3)]"<<endl
	<<"[-m <mark parts differently (default: false)]"<<endl

	<<"[-b <feature space bits size> (default: 30)]"<<endl
	<<"[-R <max radius> (default: 1)]"<<endl
	<<"[-r <min radius> (default: 0)]"<<endl
	<<"[-D <max distance relations> (default: 4)]"<<endl
	<<"[-d <min distance> (default: 0)]"<<endl
	<<"[-T <nspdk|anspdk|mnspdk|nspdk3d|alnspdk> (default: nspdk)]"<<endl
	<<"[-k output kernel matrix in output.mtx file (default: unset -> do not output)]"<<endl
	<<"[-t <hard | soft | hard_soft | mixed> as neighbourhood matching use HARD=exact matching, use SOFT=attribute matching with root identifier as radius 0 neighborhood, use HARD-SOFT=attribute matching with root identifier as full neighbourhood encoding (default: SOFT)]"<<endl
	<<"[-mink flag to set min kernel mode (default: unset->use standard dot product)]"<<endl
	<<"[-nn flag to non-normalize (default: normalization is active)]"<<endl
	<<"[-sd flag to save dot files of each RNA (default:do not save)]"<<endl
	<<"[-debug debug flag (default: false)]"<<endl;

    exit(0);
  }

  void Init(int argc, char** argv){
    vector<string> options;
    for (int i=1;i<argc; ++i ) options.push_back(argv[i]);
    for (vector<string>::iterator it=options.begin();it!=options.end();++it)
      {
	if ((*it)=="-h" || (*it)=="--help") Usage(argv[0]); 
	else if ((*it)=="-f") mInputFileName=(*(++it)); 
	else if ((*it)=="-knn") mNumNearestNeighbors=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-a") mAlpha=stream_cast<double>(*(++it)); 
	else if ((*it)=="-zl") mZThresholdLink=stream_cast<double>(*(++it)); 
	else if ((*it)=="-zc") mZThresholdCentrality=stream_cast<double>(*(++it)); 
	else if ((*it)=="-zs") mZThresholdSubtreeSize=stream_cast<double>(*(++it)); 
	else if ((*it)=="-sn") mSizeNeighborhood=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-rdl") mRootDistanceLimit=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-R") mRadiusMax=stream_cast<double>(*(++it)); 
	else if ((*it)=="-r") mRadiusMin=stream_cast<double>(*(++it)); 
	else if ((*it)=="-D") mDistanceMax=stream_cast<double>(*(++it)); 
	else if ((*it)=="-d") mDistanceMin=stream_cast<double>(*(++it)); 
	else if ((*it)=="-t") mMatchingType=(*(++it)); 
	else if ((*it)=="-v") mVerbose=true;
	else if ((*it)=="-mink") mMinKernel=true;
	else if ((*it)=="-b") mFeatureBitSize=stream_cast<int>(*(++it)); 
	else if ((*it)=="-T") mType=(*(++it)); 
	else if ((*it)=="-k") mOutputKernel=true;
	else if ((*it)=="-of") mOutputFeatures=true;
	else if ((*it)=="-nn") mNormalization=false;
	else if ((*it)=="-debug") mDebug=true;
	else if ((*it)=="-sd") mSaveDot=true;
	else if ((*it)=="-m") mMarkParts=true;
	else {cerr<<"Unrecognized parameter: "<<(*it)<<"."<<endl;throw exception();}
      }
    if (!(mMatchingType=="hard"||mMatchingType=="soft"||mMatchingType=="hard_soft"||mMatchingType=="multiview"||mMatchingType=="mixed")) {cerr<<"Wrong value for parameter: -t: "<<mMatchingType<<endl;throw exception();}
  }

public:
  string mInputFileName;
  double mRadiusMax;
  double mRadiusMin;
  double mDistanceMax;
  double mDistanceMin;
  string mMatchingType;
  bool mVerbose;
  int mFeatureBitSize;
  bool mMinKernel;
  string mType;
  bool mOutputKernel; 
  bool mNormalization;
  unsigned mNumNearestNeighbors;
  bool mSaveDot;
  double mAlpha;
  double mZThresholdLink;
  double mZThresholdCentrality;
  double mZThresholdSubtreeSize;
  unsigned mSizeNeighborhood;
  unsigned mRootDistanceLimit;
  bool mMarkParts;
  bool mOutputFeatures;
  bool mDebug;
};


class ClusterClass{
public:
  ClusterClass(NSPDK_FeatureGenerator* paFeatureGenerator, double aAlpha, bool aMarkParts):
    pmFeatureGenerator(paFeatureGenerator), 
    mAlpha(aAlpha),
    mMarkParts(aMarkParts){}

  void Input(const string& aInputFileName){
    const unsigned STEP=10;
    cerr<<"Reading data and computing features [one dot equals "<< STEP <<" graphs]..."<<endl;
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    unsigned graph_counter=0;
    while (!fin.eof()){
      GraphClass G;
      SetGraphFromFileGSPAN(fin, G);
      SVector x;
      pmFeatureGenerator->generate_feature_vector(G,x);

      //add features from stem-removed-info graph version
      SVector xs;
      GraphClass Gs=RemoveVertexLabelGivenIncidentEdgeLabel(G,"s");
      pmFeatureGenerator->generate_feature_vector(Gs,xs);
      x.add(xs);

      mDataset.push_back(x);
      graph_counter++;
      if (graph_counter%STEP==0) cerr<<".";////
      if (graph_counter%(10*STEP)==0) cerr<<"|";////
    }
    cerr<<endl;////
    fin.close();
  }

void SetGraphFromFileGSPAN(istream& in, GraphClass& oG){
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
      //status
      vector<bool> status; 
      status.push_back(true);//kernel point
      status.push_back(true);//kind
      status.push_back(true);//viewpoint
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

  GraphClass RemoveVertexLabelGivenIncidentEdgeLabel(GraphClass& aG,string aEdgeLabel){
    GraphClass G(aG);
    for (unsigned u=0;u<aG.VertexSize();u++){
      //remove label if at least one incident edge has label=aEdgeLabel
      vector<unsigned> adj=aG.GetVertexAdjacentList(u);
      for (unsigned i=0;i<adj.size();++i){
	if (aG.GetEdgeLabel(u,adj[i])==aEdgeLabel) {
	  G.SetVertexLabel(u,"X");
	  break;
	}
      }
    }
    return G;
  }

  void OutputFeatureEncoding(ostream& out){
    for (unsigned i=0;i<mDataset.size();++i)
      out<<mDataset[i];
  }

  void ComputeKernelMatrix(){
    cerr<<"Computing kernel matrix..."<<endl;
    for (unsigned i=0;i<mDataset.size();++i){
      vector<double> kernel_row;
      for (unsigned j=0;j<mDataset.size();++j){
	double kernel_ij;
	if (j<i) {
	  kernel_ij=mKernelMatrix[j][i];
	} else kernel_ij=dot(mDataset[i],mDataset[j]); 
	kernel_row.push_back(kernel_ij);
      }
      mKernelMatrix.push_back(kernel_row);
    }
  }

  void OutputKernelMatrix(ostream& out){
    for (unsigned i=0;i<mDataset.size();++i){
      for (unsigned j=0;j<mDataset.size();++j)
	out<<mKernelMatrix[i][j]<<" ";
      out<<endl;
    }
  }

  void OutputKNearestNeighbors(unsigned aNumNearestNeighbors, string aMode, ostream& out){
    unsigned effective_num_ranked=min(unsigned(mDataset.size()),aNumNearestNeighbors);
    for (unsigned i=0;i<mDataset.size();++i){
      vector<pair<double,unsigned> > similarity;
      for (unsigned j=0;j<mDataset.size();++j){
	double kernel_ij=mKernelMatrix[i][j];
	similarity.push_back(make_pair(-kernel_ij,j));//NOTE: use -kernel as sorting is done from smaller to greater
      }
      //find k most similar
      partial_sort(similarity.begin(),similarity.begin()+effective_num_ranked,similarity.end());
      //print k most similar id
      for (unsigned k=0;k<effective_num_ranked;++k)
	if (aMode=="ID") out<<similarity[k].second<<" ";
	else out<<(-similarity[k].first)<<" ";
      out<<endl;
    }
  }

  void ComputeQuickShiftClusters(double aZThresholdLink,double aZThresholdCentrality,double aZThresholdSubtreeSize, unsigned aNeigborhoodSize, unsigned aRootDistanceLimit){
    ComputeQuickShiftLinks();
    ComputeAdjacencyList();

    vector<double> subtree_size_list;
    subtree_size_list.clear();
    for (unsigned i=0;i<mDataset.size();++i)
      subtree_size_list.push_back(GetSubtreeSize(i));
    Standardize(subtree_size_list);

    cerr<<"Computing quick shift clusters..."<<endl;

    vector<double> link_strength_list;
    ofstream ofs("out_debug");//// 
    for (unsigned i=0; i<mClusterLink.size();++i){
      unsigned parent=mClusterLink[i];
      double parent_k=mKernelMatrix[i][parent];
      pair<double,double> k_stats=ComputeNeighborhoodKernelValueStatistics(i,aNeigborhoodSize);
      double avg_k=k_stats.first;
      double sd_k=k_stats.second;
      double link_strength=(parent_k-avg_k)/sd_k;
      link_strength_list.push_back(link_strength);
      ofs<<"id:"<<i<<" parent_id:"<<parent<<" parent_k:"<<parent_k<<" avg_k:"<<avg_k<<" sd_k:"<<sd_k<<" link_strength:"<<link_strength<<" centrality:"<<mCentralityList[i]<<" subtree_size:"<<GetSubtreeSize(i)<<" stdz_subtree_size:"<<subtree_size_list[i]<<" dist_from_root:"<<GetDistanceFromRoot(i)<<endl;////
    }
    Standardize(link_strength_list);

    //take as centroids those elements that:
    //1. have z score for subtree size above threshold
    //2. have z score for  above threshold
    //3. have z score for centrality above threshold
    //4. have distance from root lower than limit
    for (unsigned i=0;i<link_strength_list.size();++i){
      if (link_strength_list[i]>aZThresholdLink 
	  && mCentralityList[i]>aZThresholdCentrality 
	  && subtree_size_list[i]>aZThresholdSubtreeSize
	  && GetDistanceFromRoot(i)<=aRootDistanceLimit) {
	//remove parent link from cluster centroid (i.e. replace it with centroid id itself)
	mClusterLink[i]=i;
	//recompute tree structure
	ComputeAdjacencyList();
	subtree_size_list.clear();
	for (unsigned i=0;i<mDataset.size();++i)
	  subtree_size_list.push_back(GetSubtreeSize(i));
	Standardize(subtree_size_list);
      }
    }

    //compute cluster membership: make multimap with root element as key for each element
    for (unsigned i=0;i<mDataset.size();++i){
      mCluster.insert(make_pair(GetRoot(i),i));
    }

    //output dot file
    ofstream ofsd("out.dot");
    OutputDotFile(ofsd);
  }

  void OutputDotFile(ostream& out){
    out<<"graph{";
    for (unsigned i=0;i<mDataset.size();++i){
      double color=floor((double)i/40)/(mDataset.size()/40);
      double width=exp(mCentralityList[i]);
      out<<i<<" [shape=\"circle\" style=\"filled\" color=\""<<color<<",.4,.99\" width=\""<<width<<"\"]"<<endl;
    }
    for (unsigned i=0;i<mClusterLink.size();++i){
      unsigned j=mClusterLink[i];
      double weight=mKernelMatrix[i][j];
      out<<i<<"--"<<j<<" [weight=\""<<weight<<"\"]"<<endl;
    }
    out<<"}";
  }

  void OutputQuickShiftClusters(ostream& out){
    set<unsigned> key_list;
    for (multimap<unsigned,unsigned>::iterator it=mCluster.begin();it!=mCluster.end();++it){
      key_list.insert(it->first);
    }
    cerr<<"Determined num="<<key_list.size()<<" clusters of sizes: ";
    for (set<unsigned>::iterator it=key_list.begin();it!=key_list.end();++it){
      unsigned key=*it;
      pair<multimap<unsigned,unsigned>::iterator,multimap<unsigned,unsigned>::iterator> eqr=mCluster.equal_range(key);
      out<<"centroid: "<<key<<" cluster_members: ";
      unsigned size=0;
      for (multimap<unsigned,unsigned>::iterator jt=eqr.first;jt!=eqr.second;++jt){
	out<<jt->second<<" ";
	size++;
      }
      cerr<<size<<" ";
      out<<endl;
    }
    cerr<<endl;
  }

  void ComputeQuickShiftLinks(){
    cerr<<"Computing quick shift links..."<<endl;
    for (unsigned i=0;i<mDataset.size();++i){
      vector<pair<double,unsigned> > similarity;
      for (unsigned j=0;j<mDataset.size();++j){
	double kernel_ij=mKernelMatrix[i][j];
	similarity.push_back(make_pair(-kernel_ij,j));//NOTE: use -kernel as sorting is done from smaller to greater
      }
      //find k most similar
      sort(similarity.begin(),similarity.end());
      //find most similar id with higher centrality
      unsigned parent=i;
      for (unsigned k=0;k<similarity.size();++k){
	unsigned neighbor_id=similarity[k].second;
	if (mCentralityList[neighbor_id]>mCentralityList[i]) {parent=neighbor_id; break;}
      }
      mClusterLink.push_back(parent);
    }
    ofstream ofsd("out_complete.dot");
    OutputDotFile(ofsd);
  }

  void ComputeAdjacencyList(){
    mAdjacencyList.clear();
    for (unsigned i=0;i<mClusterLink.size();++i)
      mAdjacencyList.push_back(vector<unsigned>());
    for (unsigned i=0;i<mClusterLink.size();++i){
      unsigned child=i;
      unsigned parent=mClusterLink[i];
      if (parent!=child) mAdjacencyList[parent].push_back(child);
    }      
  }

  unsigned GetSubtreeSize(unsigned aID){
    unsigned size=1;
    if (mAdjacencyList[aID].size()!=0){
      size=1;
      for (unsigned i=0;i<mAdjacencyList[aID].size();++i)
	size+=GetSubtreeSize(mAdjacencyList[aID][i]);
    }
    return size;
  }

  unsigned GetRoot(unsigned aID){
    unsigned child=aID;
    unsigned parent=mClusterLink[child];
    while(child!=parent){
      child=parent;
      parent=mClusterLink[child];
    } 
    return parent;
  }

  unsigned GetDistanceFromRoot(unsigned aID){
    unsigned d=0;
    unsigned child=aID;
    unsigned parent=mClusterLink[child];
    while(child!=parent){
      d++;
      child=parent;
      parent=mClusterLink[child];
    } 
    return d;
  }

  pair<double,double> ComputeNeighborhoodKernelValueStatistics(unsigned aID, unsigned aNeigborhoodSize){
    vector<pair<double,unsigned> > similarity;
    for (unsigned j=0;j<mDataset.size();++j){
      double kernel_ij=mKernelMatrix[aID][j];
      similarity.push_back(make_pair(-kernel_ij,j));//NOTE: use -kernel as sorting is done from smaller to greater
    }
    partial_sort(similarity.begin(),similarity.begin()+aNeigborhoodSize,similarity.end());
    double k_avg=0;
    for (unsigned i=0;i<aNeigborhoodSize;++i){
      unsigned neighbor=similarity[i].second;
      k_avg+=mKernelMatrix[aID][neighbor];
    }
    k_avg=k_avg/aNeigborhoodSize;
    double k_sd=0;
    for (unsigned i=0;i<aNeigborhoodSize;++i){
      unsigned neighbor=similarity[i].second;
      k_sd+=(k_avg-mKernelMatrix[aID][neighbor])*(k_avg-mKernelMatrix[aID][neighbor]);
    }
    k_sd=sqrt(k_sd/aNeigborhoodSize);
    return make_pair(k_avg,k_sd);
  }

  pair<double,double> ComputeStatistics(vector<double>& aList){
    double avg=0;
    for (unsigned i=0;i<aList.size();++i){
      avg+=aList[i];
    }
    avg=avg/aList.size();
    double sd=0;
    for (unsigned i=0;i<aList.size();++i){
      sd+=(avg-aList[i])*(avg-aList[i]);
    }
    sd=sqrt(sd/aList.size());
    return make_pair(avg,sd);
  }

  void Standardize(vector<double>& aList){
    pair<double,double> stats=ComputeStatistics(aList);
    double mean=stats.first;
    double sd=stats.second;
    for (unsigned i=0;i<aList.size();++i){
      aList[i]=(aList[i]-mean)/sd;
    }
  }

  void ComputeCentrality(){
    mCentralityList.clear();
    for (unsigned i=0;i<mDataset.size();++i){
      double centrality=0;
      for (unsigned j=0;j<mDataset.size();++j)
	centrality+=mKernelMatrix[i][j];
      mCentralityList.push_back(centrality);
    }
    Standardize(mCentralityList);
  }

  void OutputCentrality(ostream& out){
    vector<pair<double,unsigned> > sorted_centrality;
    for (unsigned i=0;i<mCentralityList.size(); ++i )
      sorted_centrality.push_back(make_pair(mCentralityList[i],i));
    sort(sorted_centrality.begin(),sorted_centrality.end());
    for (unsigned i=0;i<sorted_centrality.size(); ++i )
      out<<sorted_centrality[i].second<<" "<<sorted_centrality[i].first<<endl;
  }

  void ComputeClusterRepresentatives(){
    cerr<<"Computing cluster representatives..."<<endl;
    //find most central element
    unsigned max_id=0;
    double max=mCentralityList[max_id];
    for (unsigned i=0;i<mCentralityList.size(); ++i )
      if (max<mCentralityList[i]){max=mCentralityList[i]; max_id=i;}
    

    //find highest quality element w.r.t. most central element
    unsigned max_q1_id=0;
    double max_q1=Quality(max_q1_id,max_id);
    for (unsigned i=0;i<mDataset.size(); ++i ){
      double q=Quality(i,max_id);
      if (max_q1<q){max_q1=q;max_q1_id=i;}
    }

    //find highest quality element w.r.t. previous element    
    unsigned max_q2_id=0;
    double max_q2=Quality(max_q2_id,max_q1_id);
    for (unsigned i=0;i<mDataset.size(); ++i ){
      double q=Quality(i,max_q1_id);
      if (max_q2<q){max_q2=q;max_q2_id=i;}
    }
      
    //initialize working set with previous two elements
    mCentralityList.clear();
    mClusterCentroidList.push_back(max_q1_id);
    mClusterCentroidList.push_back(max_q2_id);

    set<unsigned> working_set;
    working_set.insert(max_q1_id);
    working_set.insert(max_q2_id);
    set<unsigned> reminder_set;
    for (unsigned c=0;c<mDataset.size();++c)
      if (c!=max_q1_id && c!=max_q2_id) reminder_set.insert(c);

    //add element that has highest quality to working set (greedy procedure)
    for (unsigned r=2;r<mDataset.size();++r){//repeat n-2 times i.e. for each remaining element of the dataset
      unsigned max_qi_id=*reminder_set.begin();
      double max_qi=Quality(max_qi_id,working_set);
      for (set<unsigned>::iterator jt=reminder_set.begin(); jt!=reminder_set.end();++jt){
	unsigned q_id=*jt;
	double q=Quality(q_id,working_set);
	if (max_qi<q){max_qi=q;max_qi_id=q_id;}
      }
      mClusterCentroidList.push_back(max_qi_id);
      working_set.insert(max_qi_id);
      reminder_set.erase(max_qi_id);
    }
  }

  void OutputClusterRepresentatives(ostream& out){
    for (unsigned i=0;i<mClusterCentroidList.size();++i)
      out<<mClusterCentroidList[i]<<endl;
  }

  double Quality(unsigned i, unsigned j){
    return (mAlpha-1)*mKernelMatrix[i][j]+mAlpha*mCentralityList[i];
  }

  double Quality(unsigned i, set<unsigned>& aWorkingSet){
    double q=0;
    for (set<unsigned>::iterator jt=aWorkingSet.begin(); jt!=aWorkingSet.end();++jt){
      unsigned j=*jt;
      q+=Quality(i,j);
    }
    return q;
  }

protected:
  NSPDK_FeatureGenerator* pmFeatureGenerator;
  double mAlpha;
  unsigned mNumParts;
  unsigned mSkipNumParts;
  bool mMarkParts;
  vector<SVector> mDataset;  
  vector<double> mCentralityList;
  vector<vector<double> > mKernelMatrix;
  vector<unsigned> mClusterCentroidList;
  vector<unsigned> mClusterLink;
  multimap<unsigned,unsigned> mCluster;
  vector<vector<unsigned> > mAdjacencyList;
  vector<double> mSubtreeSizeList;
};

//---------------------------------------------------------------------------------
int main(int argc, char** argv){
  srand(time(NULL));
  LOGF<<"--------------------------------------------------------------------------------"<<endl;
  LOGF<<CREDITS<<endl;
  time_t rawtime_start; time ( &rawtime_start ); LOGF<<"Start logging: "<<asctime ( localtime ( &rawtime_start ) )<<endl;
  LOGF<<"Command line: ";for(int i=0;i<argc; ++i ) LOGF<<stream_cast<string>(argv[i])<<" ";LOGF<<endl;
  try
  {
    ParameterWrapperClass P;
    P.Init(argc,argv);
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
    else if (P.mType=="nspdk3d") pfg=&fg3d;
    else throw range_error("Unknown feature generator type:"+P.mType);

    pfg->set_flag("radius",stream_cast<string>(P.mRadiusMax));
    pfg->set_flag("distance",stream_cast<string>(P.mDistanceMax));
    pfg->set_flag("match_type",stream_cast<string>(P.mMatchingType));
    pfg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
    pfg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
    pfg->set_flag("verbosity",stream_cast<string>(P.mDebug));
    if (P.mMinKernel) pfg->set_flag("min_kernel","true");
    if (!P.mNormalization) pfg->set_flag("normalization","false");

    string ofname;

    ClusterClass C(pfg,P.mAlpha,P.mMarkParts);
    C.Input(P.mInputFileName);

    if (P.mOutputFeatures) {
      ofname=P.mInputFileName+".feature";
      ofstream ofs_f(ofname.c_str());
      C.OutputFeatureEncoding(ofs_f);      
    } else {
      C.ComputeKernelMatrix();
      ofname=P.mInputFileName+".kernel";
      ofstream ofs_k(ofname.c_str());
      C.OutputKernelMatrix(ofs_k);

      ofname=P.mInputFileName+".knn_id";
      ofstream ofs_knnid(ofname.c_str());
      C.OutputKNearestNeighbors(P.mNumNearestNeighbors, "ID", ofs_knnid);
      ofname=P.mInputFileName+".knn_val";
      ofstream ofs_knnv(ofname.c_str());
      C.OutputKNearestNeighbors(P.mNumNearestNeighbors, "VAL", ofs_knnv);

      C.ComputeCentrality();
      ofname=P.mInputFileName+".centrality";
      ofstream ofs_c(ofname.c_str());
      C.OutputCentrality(ofs_c);
      
      C.ComputeClusterRepresentatives();
      ofname=P.mInputFileName+".cluster_centroids";
      ofstream ofs_cc(ofname.c_str());
      C.OutputClusterRepresentatives(ofs_cc);


      C.ComputeQuickShiftClusters(P.mZThresholdLink, P.mZThresholdCentrality, P.mZThresholdSubtreeSize, P.mSizeNeighborhood,P.mRootDistanceLimit);
      ofname=P.mInputFileName+".quickshift_clusters";
      ofstream ofs_qsc(ofname.c_str());
      C.OutputQuickShiftClusters(ofs_qsc);
    }
  }
  catch(exception& e)
  {
    cerr<<e.what();
    LOGF<<e.what()<<endl;
  }
  time_t rawtime_end; time ( &rawtime_end );LOGF<<"End logging: "<<asctime ( localtime ( &rawtime_end ) )<<endl;
  time_t time_diff=rawtime_end-rawtime_start;LOGF<<"Time elapsed: "<<asctime ( localtime ( &time_diff ) )<<endl;
  return 0;
}


} // namespace
