
#include <iostream>

#include "sgm/Graph_boost.hh"
#include "sgm/MR_stream.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/MC_Edge.hh"
#include "sgm/Pattern.hh"

#include "dataSGMgraphs_1.icc"
#include "dataSGMgraphs_2.icc"
#include "dataSGMgraphs_3.icc"

#include "utilPrintGraph_Interface.icc"


	typedef sgm::Graph_boost<	MyGraph
								, boost::vertex_name_t
								, boost::edge_name_t
								, boost::vertex_index_t
							> GB;

void test_MC_NoEdge () {
	
	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n" 
				<<"==============================================\n" 
				<<"               sgm::MC_NoEdge  \n"
				<<"==============================================\n" 
				<<std::endl;
	
	  // load target graph
	std::string graphString;
	MyGraph targetGraph = getPattern_1( graphString );

	std::cout <<"\n == TARGET ==\n" <<std::endl;
	std::cout <<graphString <<std::endl;
	std::cout <<GB(targetGraph) <<std::endl;

	  // the target graph to fill
	MyGraph patternGraph;
	const std::string wildcard = "*";

	{
		MyGraph::vertex_descriptor v;
		MyGraph::edge_descriptor e;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), patternGraph );
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type
			edgeLabel = boost::get( boost::edge_name_t(), patternGraph );


		// node 0
		v = boost::add_vertex(patternGraph);
		nodeLabel[v] = wildcard;
		// node 1
		v = boost::add_vertex(patternGraph);
		nodeLabel[v] = wildcard;
		// node 2
		v = boost::add_vertex(patternGraph);
		nodeLabel[v] = wildcard;

		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,patternGraph), boost::vertex(1,patternGraph), patternGraph).first;
		edgeLabel[e] = wildcard;
	}

	std::cout <<"\n == PATTERN ==}\n" <<std::endl;
	GB pattern(patternGraph);
	std::cout <<pattern <<std::endl;

	sgm::SGM_vf2 sgm;
	{
		std::cout <<"\n########### matches WITHOUT CONSTRAINTS and with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}

	{
		std::cout <<"\n########### CONSTRAINED matches with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // create constraint
		sgm::MC_NoEdge noEdge = sgm::MC_NoEdge(1,2);
		std::cout <<"\n Constraint : MC_NoEdge(1,2)\n" <<std::endl;

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		constraints.push_back(&noEdge);
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n"
				<<"==============================================\n"
				<<std::endl;

}



void test_MC_EdgeLabel ( const bool allowedLabels ) {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<"==============================================\n"
				<<"               sgm::MC_EdgeLabel  \n"
				<<"==============================================\n"
				<<std::endl;

	  // load target graph
	std::string graphString;
	MyGraph targetGraph = getTarget_1( graphString );

	std::cout <<"\n == TARGET ==\n" <<std::endl;
	std::cout <<graphString <<std::endl;
	std::cout <<GB(targetGraph) <<std::endl;

	  // the target graph to fill
	MyGraph patternGraph;
	const std::string wildcard = "*";

	{
		MyGraph::vertex_descriptor v;
		MyGraph::edge_descriptor e;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), patternGraph );
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type
			edgeLabel = boost::get( boost::edge_name_t(), patternGraph );


		// node 0
		v = boost::add_vertex(patternGraph);
		nodeLabel[v] = wildcard;
		// node 1
		v = boost::add_vertex(patternGraph);
		nodeLabel[v] = wildcard;

		// edge 0 - 1
		e = boost::add_edge(boost::vertex(0,patternGraph), boost::vertex(1,patternGraph), patternGraph).first;
		edgeLabel[e] = wildcard;
	}

	std::cout <<"\n == PATTERN ==}\n" <<std::endl;
	GB pattern(patternGraph);
	std::cout <<pattern <<std::endl;

	sgm::SGM_vf2 sgm;
	{
		std::cout <<"\n########### matches WITHOUT CONSTRAINTS and with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}

	{
		std::cout <<"\n########### CONSTRAINED matches with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // create constraint
		sgm::MC_EdgeLabel::LabelSet labels;
		sgm::MC_EdgeLabel constraint = sgm::MC_EdgeLabel(0,1,labels);
		if (!allowedLabels)
			constraint = sgm::MC_EdgeLabel(0,1,labels, sgm::MC_EdgeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_EdgeLabel(0,1,{";
		for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
			std::cout <<" '"<<*l<<"'";
		if (allowedLabels)
			std::cout<<"})\n" <<std::endl;
		else
			std::cout<<"}, FORBIDDEN)\n" <<std::endl;


		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		constraints.push_back(&constraint);
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}

	{
		std::cout <<"\n########### CONSTRAINED matches with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // create constraint
		sgm::MC_EdgeLabel::LabelSet labels;
		labels.insert(wildcard);
		sgm::MC_EdgeLabel constraint = sgm::MC_EdgeLabel(0,1,labels);
		if (!allowedLabels)
			constraint = sgm::MC_EdgeLabel(0,1,labels, sgm::MC_EdgeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_EdgeLabel(0,1,{";
		for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
			std::cout <<" '"<<*l<<"'";
		if (allowedLabels)
			std::cout<<"})\n" <<std::endl;
		else
			std::cout<<"}, FORBIDDEN)\n" <<std::endl;

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		constraints.push_back(&constraint);
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}

	{
		std::cout <<"\n########### CONSTRAINED matches with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // create constraint
		sgm::MC_EdgeLabel::LabelSet labels;
		labels.insert(wildcard);
		labels.insert("-2-");
		sgm::MC_EdgeLabel constraint = sgm::MC_EdgeLabel(0,1,labels);
		if (!allowedLabels)
			constraint = sgm::MC_EdgeLabel(0,1,labels, sgm::MC_EdgeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_EdgeLabel(0,1,{";
		for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
			std::cout <<" '"<<*l<<"'";
		if (allowedLabels)
			std::cout<<"})\n" <<std::endl;
		else
			std::cout<<"}, FORBIDDEN)\n" <<std::endl;

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		constraints.push_back(&constraint);
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}

	{
		std::cout <<"\n########### CONSTRAINED matches with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // create constraint
		sgm::MC_EdgeLabel::LabelSet labels;
		labels.insert("-2-");
		sgm::MC_EdgeLabel constraint = sgm::MC_EdgeLabel(0,1,labels);
		if (!allowedLabels)
			constraint = sgm::MC_EdgeLabel(0,1,labels, sgm::MC_EdgeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_EdgeLabel(0,1,{";
		for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
			std::cout <<" '"<<*l<<"'";
		if (allowedLabels)
			std::cout<<"})\n" <<std::endl;
		else
			std::cout<<"}, FORBIDDEN)\n" <<std::endl;

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		constraints.push_back(&constraint);
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}

	{
		std::cout <<"\n########### CONSTRAINED matches with wildcard '" <<wildcard <<"' ###########\n" <<std::endl;
		sgm::MR_stream mr(std::cout);

		  // create constraint
		sgm::MC_EdgeLabel::LabelSet labels;
		labels.insert("-1-");
		labels.insert("-2-");
		sgm::MC_EdgeLabel constraint = sgm::MC_EdgeLabel(0,1,labels);
		if (!allowedLabels)
			constraint = sgm::MC_EdgeLabel(0,1,labels, sgm::MC_EdgeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_EdgeLabel(0,1,{";
		for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
			std::cout <<" '"<<*l<<"'";
		if (allowedLabels)
			std::cout<<"})\n" <<std::endl;
		else
			std::cout<<"}, FORBIDDEN)\n" <<std::endl;

		  // setup pattern without constraint
		sgm::Pattern_Interface::ConstraintVec constraints;
		constraints.push_back(&constraint);
		sgm::Pattern toFind = sgm::Pattern(pattern, constraints, wildcard);
		  // perform matching
		sgm.findMatches( toFind, GB(targetGraph), mr, UINT_MAX );
	}


	std::cout	<<"\n"
				<<"===============  TEST END  ===================\n"
				<<"==============================================\n"
				<<std::endl;

}

int main(int argc, char** argv) {

	  // test MC_NoEdge
	test_MC_NoEdge();

	  // test MC_EdgeLabel
	test_MC_EdgeLabel( true );
	test_MC_EdgeLabel( false );

	
	return 0;
}


