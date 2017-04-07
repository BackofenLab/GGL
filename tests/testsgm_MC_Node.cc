
#include <iostream>

#include "sgm/Graph_boost.hh"
#include "sgm/MR_stream.hh"
#include "sgm/SGM_vf2.hh"
#include "sgm/MC_Node.hh"
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



void test_MC_NodeLabel ( const bool allowedLabels ) {

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<"==============================================\n"
				<<"               sgm::MC_NodeLabel  \n"
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
		nodeLabel[v] = "C";

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
		sgm::MC_NodeLabel::LabelSet labels;
		sgm::MC_NodeLabel constraint = sgm::MC_NodeLabel(0,labels);
		if (!allowedLabels)
			constraint = sgm::MC_NodeLabel(0,labels, sgm::MC_NodeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_NodeLabel(0,{";
		for (sgm::MC_NodeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
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
		sgm::MC_NodeLabel::LabelSet labels;
		labels.insert(wildcard);
		sgm::MC_NodeLabel constraint = sgm::MC_NodeLabel(0,labels);
		if (!allowedLabels)
			constraint = sgm::MC_NodeLabel(0,labels, sgm::MC_NodeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_NodeLabel(0,{";
		for (sgm::MC_NodeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
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
		sgm::MC_NodeLabel::LabelSet labels;
		labels.insert(wildcard);
		labels.insert("N");
		sgm::MC_NodeLabel constraint = sgm::MC_NodeLabel(0,labels);
		if (!allowedLabels)
			constraint = sgm::MC_NodeLabel(0,labels, sgm::MC_NodeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_NodeLabel(0,{";
		for (sgm::MC_NodeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
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
		sgm::MC_NodeLabel::LabelSet labels;
		labels.insert("N");
		sgm::MC_NodeLabel constraint = sgm::MC_NodeLabel(0,labels);
		if (!allowedLabels)
			constraint = sgm::MC_NodeLabel(0,labels, sgm::MC_NodeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_NodeLabel(0,{";
		for (sgm::MC_NodeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
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
		sgm::MC_NodeLabel::LabelSet labels;
		labels.insert("O");
		labels.insert("N");
		sgm::MC_NodeLabel constraint = sgm::MC_NodeLabel(0,labels);
		if (!allowedLabels)
			constraint = sgm::MC_NodeLabel(0,labels, sgm::MC_NodeLabel::FORBIDDEN);
		std::cout <<"\n Constraint : MC_NodeLabel(0,{";
		for (sgm::MC_NodeLabel::LabelSet::const_iterator l=labels.begin(); l!=labels.end();++l)
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

void test_MC_NodeAdjacency( const MyGraph & targetGraph, const sgm::MC_NodeAdjacency::MC_Operator op, const size_t count ) {

	std::string opString = "";
	switch(op) {
	case sgm::MC_NodeAdjacency::MC_EQ : opString = "="; break;
	case sgm::MC_NodeAdjacency::MC_L : opString = "<"; break;
	case sgm::MC_NodeAdjacency::MC_LQ : opString = "<="; break;
	case sgm::MC_NodeAdjacency::MC_G : opString = ">"; break;
	case sgm::MC_NodeAdjacency::MC_GQ : opString = ">="; break;
	default: std::cout <<" unknown operator "<<op<<std::endl;
	}

	std::cout	<<"\n"
				<<"==============  BEGIN TEST  ==================\n"
				<<"==============================================\n"
				<<"               sgm::MC_NodeLabel  \n"
				<<"==============================================\n"
				<<std::endl;

	  // print target graph
	std::cout <<"\n == TARGET ==\n" <<std::endl;
	std::cout <<GB(targetGraph) <<std::endl;

	  // the target graph to fill
	MyGraph patternGraph;
	const std::string wildcard = "*";

	{
		MyGraph::vertex_descriptor v;
		MyGraph::edge_descriptor e;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), patternGraph );


		// node 0
		v = boost::add_vertex(patternGraph);
		nodeLabel[v] = wildcard;
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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert(wildcard);
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert(wildcard);
		nodeLabels.insert("C");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		nodeLabels.insert("N");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		edgeLabels.insert("-2-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert(wildcard);
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		nodeLabels.insert("N");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert(wildcard);
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert(wildcard);
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert(wildcard);
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		edgeLabels.insert("-2-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		edgeLabels.insert("-2-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		nodeLabels.insert("N");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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
		sgm::MC_NodeAdjacency::LabelSet nodeLabels;
		nodeLabels.insert("C");
		nodeLabels.insert("N");
		sgm::MC_NodeAdjacency::LabelSet edgeLabels;
		edgeLabels.insert("-1-");
		edgeLabels.insert("-2-");
		std::cout <<"\n Constraint : MC_NodeAdjacency(0, "<<opString<<", " <<count<<", {";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=nodeLabels.begin(); l!=nodeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout <<" },{";
		for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=edgeLabels.begin(); l!=edgeLabels.end();++l)
			std::cout <<" '"<<*l<<"'";
		std::cout<<" })\n" <<std::endl;
		sgm::MC_NodeAdjacency constraint = sgm::MC_NodeAdjacency(0,op,count,nodeLabels,edgeLabels);


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


	  // test MC_NodeLabel
	test_MC_NodeLabel( true );
	test_MC_NodeLabel( false );


	  // test MC_NodeAdjacency
	std::string targetString;
	MyGraph targetGraph = getTarget_1( targetString );

	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 0 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 1 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 2 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 3 );
	
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_L, 0 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_L, 1 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_L, 2 );

	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_G, 1 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_G, 2 );

	// TODO test loop and parallel edge cases
	targetGraph = MyGraph();
	{
		MyGraph::vertex_descriptor v;
		MyGraph::edge_descriptor e;

		boost::property_map< MyGraph, boost::vertex_name_t >::type
			nodeLabel = boost::get( boost::vertex_name_t(), targetGraph );
		boost::property_map< MyGraph, boost::vertex_index_t >::type
			nodeIndex = boost::get( boost::vertex_index_t(), targetGraph );
		  // get level property class
		boost::property_map< MyGraph, boost::edge_name_t >::type
			edgeLabel = boost::get( boost::edge_name_t(), targetGraph );

		size_t from = 0, to = 0;

		// node 0
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		from = nodeIndex[v];
		to = nodeIndex[v];
		// edge 0 - 0
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";

		// node 1
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		from = nodeIndex[v];
		to = nodeIndex[v];
		// edge 1 - 1
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";

		// node 2
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		from = nodeIndex[v];
		to = nodeIndex[v];
		// edge 2 - 2
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-2-";
		// node 3
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		to = nodeIndex[v];
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";

		// node 2
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		from = nodeIndex[v];
		to = nodeIndex[v];
		// edge 2 - 2
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-2-";
		// node 3
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "N";
		to = nodeIndex[v];
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";

		// node 2
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		from = nodeIndex[v];
		to = nodeIndex[v];
		// edge 2 - 2
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-2-";
		// node 3
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		to = nodeIndex[v];
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-2-";

		// node 2
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "C";
		from = nodeIndex[v];
		to = nodeIndex[v];
		// edge 2 - 2
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-2-";
		// node 3
		v = boost::add_vertex(targetGraph);
		nodeLabel[v] = "N";
		to = nodeIndex[v];
		// edge 2 - 3
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-1-";
		e = boost::add_edge(boost::vertex(from,targetGraph), boost::vertex(to,targetGraph), targetGraph).first;
		edgeLabel[e] = "-2-";


	}

	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 0 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 1 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 2 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 3 );
	test_MC_NodeAdjacency( targetGraph, sgm::MC_NodeAdjacency::MC_EQ, 4 );

	return 0;
}


