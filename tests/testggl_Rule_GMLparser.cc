
#include "ggl/Rule_GMLparser.hh"

#include "sgm/Graph_boost.hh"


#include "utilPrintRule.icc"


std::string gmlRule1 = std::string(
"rule [\n\
		ruleID \"gmlRule1\"\
        context [\n\
                node [ id 0 label \"C\" ]\n\
                node [ id 1 label \"C\" ]\n\
                node [ id 2 label \"C\" ]\n\
        ]\n\
        left [\n\
                edge [ source 1 target 2 label \"-\" ]\n\
        ]\n\
        right [\n\
                edge [ source 0 target 1 label \"-\" ]\n\
                edge [ source 2 target 0 label \"-\" ]\n\
        ]\n\
]");

const std::string Rule1 = std::string("\
  3(D) -1- 1(B)               3(E) -1- 1(B) \n\
   |      /                    |        |   \n\
   2    3           ==>        4        3   \n\
   |  /                        |        |   \n\
  2(C) -0- 0(A)               2(C)     4(D) \n\
");

const std::string Rule1gml = std::string(
"rule [\n\
		ruleID \"Rule1gml\"\
        context [  \n\
                node [ id 1 label \"B\" ]  \n\
                node [ id 2 label \"C\" ]  \n\
                edge [ source 1 target 3 label \"-1-\" ]  \n\
        ]  \n\
        left [  \n\
                node [ id 0 label \"A\" ]  \n\
                node [ id 3 label \"D\" ]  \n\
                edge [ source 0 target 2 label \"-0-\" ]  \n\
                edge [ source 1 target 2 label \"-3-\" ]  \n\
                edge [ source 2 target 3 label \"-2-\" ]  \n\
        ]  \n\
        right [  \n\
                node [ id 3 label \"E\" ]\n\
                node [ id 4 label \"D\" ]\n\
                edge [ source 2 target 3 label \"-4-\" ]  \n\
                edge [ source 1 target 4 label \"-3-\" ]  \n\
        ]  \n\
]\
");



const std::string CompactRule1gml = std::string(
"graph [\n\
	node [ id 0 label \"A|\" ]  \n\
	node [ id 1 label \"B\" ]  \n\
	node [ id 2 label \"C\" ]  \n\
	node [ id 3 label \"D|E\" ]\n\
	node [ id 4 label \"|D\" ]\n\
	edge [ source 0 target 2 label \"-0-|\" ]  \n\
	edge [ source 1 target 3 label \"-1-\" ]  \n\
	edge [ source 1 target 2 label \"-3-|\" ]  \n\
	edge [ source 2 target 3 label \"-2-|-4-\" ]  \n\
	edge [ source 1 target 4 label \"|-3-\" ]  \n\
]\
");


const std::string Rule2 = std::string("\
  0(C)   1(C) * 2(C)               0(C) - 1(C) - 2(C) \n\
                                                      \n\
   ||                     ==>       |             |   \n\
                                                      \n\
  5(*) - 4(C) = 3(C)               5(*) = 4(C) - 3(C) \n\
");

const std::string Rule2gml = std::string("\
rule [\n\
 ruleID \"Rule2gml\"\
 context [\n\
   node [ id 0 label \"C\" ]\n\
   node [ id 1 label \"C\" ]\n\
   node [ id 2 label \"C\" ]\n\
   node [ id 3 label \"C\" ]\n\
   node [ id 4 label \"C\" ]\n\
   node [ id 5 label \"*\" ]\n\
 ]\n\
 left [\n\
   edge [ source 1 target 2 label \"*\" ]\n\
   edge [ source 3 target 4 label \"=\" ]\n\
   edge [ source 4 target 5 label \"-\" ]\n\
   edge [ source 5 target 0 label \"=\" ]\n\
 ]\n\
 right [\n\
   edge [ source 0 target 1 label \"-\" ]\n\
   edge [ source 1 target 2 label \"-\" ]\n\
   edge [ source 2 target 3 label \"-\" ]\n\
   edge [ source 3 target 4 label \"-\" ]\n\
   edge [ source 4 target 5 label \"=\" ]\n\
   edge [ source 5 target 0 label \"-\" ]\n\
 ]\n\
 constrainAdj [\n\
   id 0\n\
   op =\n\
   count 5\n\
   nodeLabels [ label \"H\" label \"X\" ]\n\
   edgeLabels [ label \"-\" ]\n\
 ]\n\
 constrainAdj [\n\
   id 3\n\
   op <\n\
   count 2\n\
   edgeLabels [ label \"=\" ]\n\
 ]\n\
 constrainAdj [\n\
   id 2\n\
   op >\n\
   count 0\n\
 ]\n\
 constrainNode [\n\
   id 5\n\
   nodeLabels [ label \"H\" label \"X\" ]\n\
 ]\n\
 constrainNode [\n\
   id 5\n\
   op =\n\
   nodeLabels [ label \"H\" label \"X\" ]\n\
 ]\n\
 constrainNode [\n\
   id 5\n\
   op !\n\
   nodeLabels [ label \"Z\" ]\n\
 ]\n\
 constrainEdge [\n\
   source 1\n\
   target 2\n\
   edgeLabels [ label \"=\" label \"#\" ]\n\
 ]\n\
 constrainEdge [\n\
   source 1\n\
   target 2\n\
   op =\n\
   edgeLabels [ label \"=\" label \"#\" ]\n\
 ]\n\
 constrainEdge [\n\
   source 1\n\
   target 2\n\
   op !\n\
   edgeLabels [ label \"-\" ]\n\
 ]\n\
 constrainNoEdge [\n\
   source 0\n\
   target 2\n\
 ]\n\
]\
");


const std::string Rule3 = std::string("\
  1(A) --- 2(B)   -->  2(B)+edges(A)\n\
");

const std::string Rule3gml = std::string(
"rule [\n\
        ruleID \"Rule3gml\"\
        context [  \n\
                node [ id 2 label \"B\" ]  \n\
        ]  \n\
        left [  \n\
                node [ id 1 label \"A\" ]  \n\
                edge [ source 1 target 2 label \"---\" ]  \n\
        ]  \n"
"        wildcard \"A\"\n"
"        copyAndPaste [ source 1 id 2 ]\n"
"        copyAndPaste [ source 1 id 2 target 2 ]\n"
"        copyAndPaste [ source 1 id 2 edgeLabels [ label \"+++\" ] ]\n"
"        copyAndPaste [ source 1 id 2 target 2 edgeLabels [ label \"+++\" ] ]\n"
"]"
);

void printConstraints( const ggl::Rule rule );
void printCopyAndPaste( const ggl::Rule rule );
void printWildcard( const ggl::Rule rule );

int main( int argc, char** argv ) {
	
	std::vector<std::string> rules;
	std::vector<std::string> rulesGML;
	
	rules.push_back(Rule1);
	rulesGML.push_back(Rule1gml);
	rules.push_back(Rule2);
	rulesGML.push_back(Rule2gml);
	rules.push_back(Rule3);
	rulesGML.push_back(Rule3gml);
	
	
	for (size_t i=0; i<rules.size(); i++) {
		
		std::cout	<<"\n======= RULE GML TO PARSE =======\n" 
					<<rules[i] 
					<<"\n=================================\n"
					<<rulesGML[i]
					<<"\n================================="
					<<std::endl;
		
		try {
			std::pair<ggl::Rule, int> ret = ggl::Rule_GMLparser::parseRule( rulesGML[i] );


			if (ret.second != -1) {
				std::cout <<" Parsing error at position " <<ret.second <<std::endl;
			}

			std::cout	<<"\n========= PARSED RULE ===========" <<std::endl;
			printRule(ret.first);
			printLeftSidePattern(ret.first);
			printRightSidePattern(ret.first);
			printConstraints(ret.first);
			printCopyAndPaste(ret.first);
			printWildcard(ret.first);
			std::cout	<<"\n=================================" <<std::endl;
		} catch (ggl::Rule_GML_error &ex) {
			  // catch GML parsing errors
			std::cout	<<" EXPECTION RAISED : Rule parsing error in rule specification :\n"
						<<ex.what()
						<<std::endl;
		}

	}
	

	/////////// TEST COMPACTED RULE PARSING

	rules.clear();
	rulesGML.clear();

	rules.push_back(Rule1);
	rulesGML.push_back(CompactRule1gml);


	for (size_t i=0; i<rules.size(); i++) {

		std::cout	<<"\n== RULE COMPACTED GML TO PARSE ==\n"
					<<rules[i]
					<<"\n=================================\n"
					<<rulesGML[i]
					<<"\n================================="
					<<std::endl;

		try {
			std::pair<ggl::Rule, int> ret = ggl::Rule_GMLparser::parseCompactRule( rulesGML[i] );

			if (ret.second != -1) {
				std::cout <<" Parsing error at position " <<ret.second <<std::endl;
			}

			std::cout	<<"\n========= PARSED RULE ===========" <<std::endl;
			printRule(ret.first);
			printLeftSidePattern(ret.first);
			printRightSidePattern(ret.first);
			printConstraints(ret.first);
			printCopyAndPaste(ret.first);
			printWildcard(ret.first);
			std::cout	<<"\n=================================" <<std::endl;
		} catch (ggl::Rule_GML_error &ex) {
			  // catch GML parsing errors
			std::cout	<<" EXPECTION RAISED : Rule parsing error in rule specification :\n"
						<<ex.what()
						<<std::endl;
		}

	}

	
	return 0;
}

void printConstraints( const ggl::Rule rule )
{
	std::cout <<"\n number of constraints = " <<rule.getLeftSide().constraints.size() <<std::endl;

	typedef std::vector< sgm::Pattern_Interface::Match_Constraint* > CV;
	const CV & cv = rule.getLeftSide().constraints;
	for (CV::const_iterator c = cv.begin(); c != cv.end(); ++c ) {
		const sgm::MC_NodeLabel * cur = dynamic_cast< const sgm::MC_NodeLabel * >(*c);
		{
		if (cur != NULL) {
			std::cout <<" + MC_NodeLabel "<<std::endl;
			std::cout <<"   id " <<cur->constrainedNodeID <<std::endl;
			std::cout <<"   op " <<(cur->compareType==sgm::MC_NodeLabel::ALLOWED?"=":"!") <<std::endl;
			std::cout <<"   nodeLabels = " ;
			for (sgm::MC_NodeLabel::LabelSet::const_iterator l=cur->nodeLabels.begin(); l!=cur->nodeLabels.end(); ++l)
				std::cout <<" '"<<*l<<"'" ;
			std::cout <<std::endl;
			continue;
		}}{
		const sgm::MC_EdgeLabel * cur = dynamic_cast< const sgm::MC_EdgeLabel * >(*c);
		if (cur != NULL) {
			std::cout <<" + MC_EdgeLabel "<<std::endl;
			std::cout <<"   from " <<cur->constrainedFromID <<std::endl;
			std::cout <<"   to   " <<cur->constrainedToID <<std::endl;
			std::cout <<"   op " <<(cur->compareType==sgm::MC_EdgeLabel::ALLOWED?"=":"!") <<std::endl;
			std::cout <<"   edgeLabels = " ;
			for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=cur->edgeLabels.begin(); l!=cur->edgeLabels.end(); ++l)
				std::cout <<" '"<<*l<<"'" ;
			std::cout <<std::endl;
			continue;
		}}{
		const sgm::MC_NoEdge * cur = dynamic_cast< const sgm::MC_NoEdge * >(*c);
		if (cur != NULL) {
			std::cout <<" + MC_NoEdge "<<std::endl;
			std::cout <<"   from " <<cur->constrainedFromID <<std::endl;
			std::cout <<"   to   " <<cur->constrainedToID <<std::endl;
			continue;
		}}{
		const sgm::MC_NodeAdjacency * cur = dynamic_cast< const sgm::MC_NodeAdjacency * >(*c);
		if (cur != NULL) {
			std::cout <<" + MC_NodeAdjacency "<<std::endl;
			std::cout <<"   id " <<cur->constrainedNodeID <<std::endl;
			std::cout <<"   op " ;
			switch(cur->op) {
			case sgm::MC_NodeAdjacency::MC_EQ : std::cout <<"="; break;
			case sgm::MC_NodeAdjacency::MC_G : std::cout <<">"; break;
			case sgm::MC_NodeAdjacency::MC_GQ : std::cout <<">="; break;
			case sgm::MC_NodeAdjacency::MC_L : std::cout <<"<"; break;
			case sgm::MC_NodeAdjacency::MC_LQ : std::cout <<"<="; break;
			default : std::cout <<"unknown operator type";
			}
			std::cout <<std::endl;
			std::cout <<"   count " <<cur->count <<std::endl;
			std::cout <<"   nodeLabels = " ;
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=cur->nodeLabels.begin(); l!=cur->nodeLabels.end(); ++l)
				std::cout <<" '"<<*l<<"'" ;
			std::cout <<std::endl;
			std::cout <<"   edgeLabels = " ;
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=cur->edgeLabels.begin(); l!=cur->edgeLabels.end(); ++l)
				std::cout <<" '"<<*l<<"'" ;
			std::cout <<std::endl;

			continue;
		}}
		std::cout <<" !!! unknown constraint type !!! " <<std::endl;
	}

}

void printCopyAndPaste( const ggl::Rule rule )
{
	std::cout <<"\n number of nodes with copy-and-paste operations = " <<rule.getCopyAndPasteOperations().size() <<"\n"<<std::endl;
	for ( ggl::Rule::CopyAndPasteOperations::const_iterator n=rule.getCopyAndPasteOperations().begin(); n!=rule.getCopyAndPasteOperations().end(); n++) {
		for (ggl::Rule::CopyAndPasteOperations::mapped_type::const_iterator cnp=n->second.begin(); cnp!=n->second.end(); cnp++) {
			std::cout <<" copy-and-paste : source("<<cnp->source <<") pasteID("<<cnp->pasteID<<") target(";
			if (cnp->target != (size_t)INT_MAX)
				std::cout <<cnp->target;
			std::cout <<") edgeLabels (";
			for (ggl::Rule::RuleCnP::LabelSet::const_iterator l=cnp->edgeLabels.begin(); l!=cnp->edgeLabels.end(); ++l)
				std::cout <<" '"<<*l<<"'";
			std::cout <<")"<<std::endl;
		}
	}
}

void printWildcard( const ggl::Rule rule )
{
	std::cout <<"\n wildcard pattern = " <<(rule.getUsedWildcard()==NULL?"NULL" : "'"+*(rule.getUsedWildcard())+"'") <<std::endl;
}
