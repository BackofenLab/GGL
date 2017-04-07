
#include <iostream>
#include <algorithm>

#include "ggl/Graph.hh"

#include "dataTargetGraph_1.icc"

#include <boost/graph/breadth_first_search.hpp>

#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_utility.hpp>

template <class GRAPH>
struct cycle_reporter 
  : public boost::base_visitor<cycle_reporter<GRAPH> >
{
  typedef boost::on_gray_target event_filter;
//  typedef boost::on_black_target event_filter;
//  typedef boost::on_examine_edge event_filter;

  cycle_reporter(const std::vector<typename GRAPH::vertex_descriptor>& pred_)
   : pred(pred_)
  { }

  template <class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
    std::cout <<" found cycle closing : ("<<boost::source(e, g)<<", " << boost::target(e, g) <<std::endl;
    typename GRAPH::vertex_descriptor v = boost::target(e,g);
	std::cout <<" " <<(size_t)v;
    while( pred[v] != v ){
    	v = pred[v];
    	std::cout <<" " <<(size_t)v;
    } 
    std::cout <<" +";
    v = boost::source(e,g);
	std::cout <<" " <<(size_t)v;
    while( pred[v] != v ){
    	v = pred[v];
    	std::cout <<" " <<(size_t)v;
    } 
    std::cout <<std::endl;
  }
protected:
	const std::vector<typename GRAPH::vertex_descriptor>& pred;
};

template <class GRAPH>
struct edge_reporter 
  : public boost::base_visitor<edge_reporter<GRAPH> >
{
//  typedef boost::on_black_target event_filter;
//  typedef boost::on_examine_edge event_filter;
  typedef boost::on_tree_edge event_filter;

  edge_reporter(const std::vector<typename GRAPH::vertex_descriptor>& pred_)
   : pred(pred_)
  { }

  template <class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
    std::cout <<" next edge : ("<<boost::source(e, g)<<", " << boost::target(e, g) <<std::endl;
  }
protected:
	const std::vector<typename GRAPH::vertex_descriptor>& pred;
};

template <class GRAPH>
struct shortest_path_reporter 
  : public boost::base_visitor<shortest_path_reporter<GRAPH> >
{
  typedef boost::on_tree_edge event_filter;

  shortest_path_reporter(	
		  	const std::vector<typename GRAPH::vertex_descriptor>& pred_
		  	, const typename GRAPH::vertex_descriptor & target_
		  	, std::vector<typename GRAPH::vertex_descriptor>& path_)
   : pred(pred_), target(target_), path(path_)
  { 
  }

  template <class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
	  if (boost::target(e,g) == target) {
		    typename GRAPH::vertex_descriptor v = boost::source(e,g);
		    while( pred[v] != v ){
		    	  // add to path
		    	path.push_back(v);
		    	v = pred[v];
	    	}
	  }
  }
protected:
	const std::vector<typename GRAPH::vertex_descriptor>& pred;
	const typename GRAPH::vertex_descriptor & target;
	std::vector<typename GRAPH::vertex_descriptor>& path;
};




//inline cycle_reporter
//report_cycle() {
//  return cycle_reporter();
//}

int main(int argc, char** argv) {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"              ggl::CE_Horton  \n" 
				<<"      detection of all cycles in graphs \n" 
				<<"==============================================\n" 
				<<std::endl;
	
	std::string graphString = "";
	
	std::cout <<"-->create ggl::Graph = getTarget_1()" <<std::endl;
	
	ggl::Graph t = getTarget_1(graphString);

	
	std::cout <<"\n" <<"filled graph is :\n\n" <<graphString <<std::endl;
	
	ggl::Graph::vertex_iterator vi,v_end;
	ggl::Graph::edge_iterator ei,e_end;
	
	for (boost::tie(ei,e_end)=boost::edges(t); ei!=e_end; ei++)  {
		ggl::Graph::vertex_descriptor i = boost::source(*ei,t);
		ggl::Graph::vertex_descriptor j = boost::target(*ei,t);
		for (boost::tie(vi,v_end)=boost::vertices(t); vi!=v_end; vi++)  {
			if ( *vi != i && *vi != j ) {
				std::vector<ggl::Graph::vertex_descriptor> path_v_i;
				path_v_i.push_back(*vi);
				{				
					  // the predecessor mapping
					std::vector<ggl::Graph::vertex_descriptor> p(boost::num_vertices(t));	
					  // The source vertex maps to itself
					p[*vi] = *vi;
					  // find path from *vi --> i
					boost::breadth_first_search
					(t, *vi, 
					boost::visitor(boost::make_bfs_visitor(
							std::make_pair(
									boost::record_predecessors(&p[0],boost::on_tree_edge())
									, shortest_path_reporter<ggl::Graph>(p,i,path_v_i)))) );
				}
				path_v_i.push_back(i);
				std::vector<ggl::Graph::vertex_descriptor> path_v_j;
				{				
					  // the predecessor mapping
					std::vector<ggl::Graph::vertex_descriptor> p(boost::num_vertices(t));	
					  // The source vertex maps to itself
					p[*vi] = *vi;
					  // find path from *vi --> j
					boost::breadth_first_search
					(t, *vi, 
					boost::visitor(boost::make_bfs_visitor(
							std::make_pair(
									boost::record_predecessors(&p[0],boost::on_tree_edge())
									, shortest_path_reporter<ggl::Graph>(p,j,path_v_j)))) );
				}
				path_v_j.push_back(j);
				
				bool isUnique = true;
				for (size_t ix=0; isUnique && ix<path_v_i.size(); ix++) {
					isUnique = std::find(path_v_j.begin(), path_v_j.end(), path_v_i[ix]) == path_v_j.end();
				}
				if (isUnique) {
					  // reverse last part
					std::cout <<"  = next cycle = ";
					copy(path_v_i.begin(), path_v_i.end(), std::ostream_iterator<ggl::Graph::vertex_descriptor>(std::cout, " "));
					copy(path_v_j.rbegin(), path_v_j.rend(), std::ostream_iterator<ggl::Graph::vertex_descriptor>(std::cout, " "));
					std::cout <<std::endl;
				}
			}
		}
	}
	
	

	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n" 
				<<"==============================================\n"
				<<std::endl;
	
	return 0;
}
