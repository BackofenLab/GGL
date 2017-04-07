
#include "sgm/Pattern.hh"

#include <stack>
#include <cassert>
#include <climits>
#include <set>



namespace sgm {


////////////////////////////////////////////////////////////////////////////

	Pattern_Interface
	:: ~Pattern_Interface()
	{
	}


////////////////////////////////////////////////////////////////////////////

	  // Inequality comparison
	  // @param toCompare the Pattern to compare to
	  // @return true if both interfaces describe different patterns
	bool
	Pattern_Interface
	::operator==(const Pattern_Interface& p ) const
	{
		  // check if graph pattern is equal
		  // and same number of match constraints
		bool isEqual = this->getPatternGraph().operator==(p.getPatternGraph())
				&& this->getConstraints().size() == p.getConstraints().size();

		  // compare the constraints (NOTE: have to show the same order)
		for (size_t i=0; isEqual && i<this->getConstraints().size(); ++i)
		{
			  // pairwise constraint comparison
			isEqual = *(this->getConstraints().at(i)) == *(p.getConstraints().at(i));
		}

		return isEqual &&
				(this->getUsedWildcard() == p.getUsedWildcard()
					|| (this->getUsedWildcard() != NULL
							&& p.getUsedWildcard() != NULL
							&& *(this->getUsedWildcard()) == *(p.getUsedWildcard()) ));
	}


////////////////////////////////////////////////////////////////////////////

	  // Inequality comparison
	  // @param toCompare the Pattern to compare to
	  // @return true if both interfaces describe different patterns
	bool
	Pattern_Interface
	::operator!=(const Pattern_Interface& p ) const
	{
		return !(this->operator==(p));
	}


////////////////////////////////////////////////////////////////////////////


} // namespace sgm


namespace sgm {


////////////////////////////////////////////////////////////////////////////

	// constructs a pattern without additional constraints
	// @param graph the graph to represent
	Pattern
	::Pattern( const Graph_Interface& graph_ )
	 :	graph(&graph_)
		, matchConstraints()
		, usedWildcard(NULL)
	{
	}


////////////////////////////////////////////////////////////////////////////

	// constructs a pattern without additional constraints
	// @param graph the graph to represent
	Pattern
	::Pattern( const Graph_Interface& graph_,
				const std::string & wildcardToUse)
	 :	graph(&graph_)
		, matchConstraints()
		, usedWildcard(new std::string(wildcardToUse))
	{
	}


////////////////////////////////////////////////////////////////////////////

	// constructs a pattern with additional constraints
	// @param graph the graph to represent
	// @param matchConstraints the additional constraints to be fulfilled
	Pattern
	::Pattern( const Graph_Interface& graph_
			, const Pattern_Interface::ConstraintVec & matchConstraints_ )
	 :	graph(&graph_)
		, matchConstraints( copyConstraintVec(matchConstraints_) )
		, usedWildcard(NULL)
	{

	}


////////////////////////////////////////////////////////////////////////////

	// constructs a pattern with additional constraints
	// @param graph the graph to represent
	// @param matchConstraints the additional constraints to be fulfilled
	Pattern
	::Pattern( const Graph_Interface& graph_
			, const Pattern_Interface::ConstraintVec & matchConstraints_
			, const std::string & wildcardToUse )
	 :	graph(&graph_)
		, matchConstraints( copyConstraintVec(matchConstraints_) )
		, usedWildcard(new std::string(wildcardToUse))
	{

	}


////////////////////////////////////////////////////////////////////////////

	// copy construction
	// @param the pattern to make this object a copy of
	Pattern
	::Pattern( const Pattern& toCopy )
	 :	graph( toCopy.graph )
		, matchConstraints( copyConstraintVec(toCopy.matchConstraints) )
		, usedWildcard( toCopy.usedWildcard==NULL ? NULL : new std::string(*usedWildcard))
	{
	}


////////////////////////////////////////////////////////////////////////////

	Pattern
	:: ~Pattern()
	{
		  // remove local MatchConstraint instances
		for (size_t i=0; i<matchConstraints.size(); ++i) {
			delete(matchConstraints[i]);
		}
		matchConstraints.clear();
		  // remove used wildcard object
		if (usedWildcard != NULL)
			delete usedWildcard;
	}


////////////////////////////////////////////////////////////////////////////


	const Graph_Interface&
	Pattern
	::getPatternGraph( void ) const
	{
		  // access to the pattern graph
		return *graph;
	}



////////////////////////////////////////////////////////////////////////////


	const Pattern_Interface::ConstraintVec&
	Pattern
	::getConstraints( void ) const
	{
		  // access to matching constraints
		return matchConstraints;
	}


////////////////////////////////////////////////////////////////////////////


	const std::string*
	Pattern
	::getUsedWildcard( void ) const
	{
		  // access to the wildcard to be used for matching
		return usedWildcard;
	}


////////////////////////////////////////////////////////////////////////////


	Pattern_Interface::ConstraintVec
	Pattern
	::copyConstraintVec( const ConstraintVec & toCopy )
	{
		  // create a deep copy of toCopy, i.e. clone each referenced object
		ConstraintVec copy( toCopy.size(), NULL );
		for (size_t i=0; i<copy.size(); ++i)  {
			copy[i] = toCopy.at(i)->clone();
		}
		return copy;
	}


////////////////////////////////////////////////////////////////////////////


	Pattern&
	Pattern
	::operator=(const sgm::Pattern& toCopy )
	{
		  // remove local MatchConstraint instances
		for (size_t i=0; i<matchConstraints.size(); ++i) {
			delete(matchConstraints[i]);
		}
		matchConstraints.clear();
		  // remove used wildcard object
		if (usedWildcard != NULL)
			delete usedWildcard;

		 // copy data
		graph = toCopy.graph;
		matchConstraints = copyConstraintVec(toCopy.matchConstraints);
		usedWildcard = toCopy.usedWildcard==NULL ? NULL : new std::string(*toCopy.usedWildcard);

		return *this;
	}


////////////////////////////////////////////////////////////////////////////



} // namespace sgm

std::ostream&
operator <<( std::ostream & out, const sgm::Pattern& g )
{
		out <<" pattern graph :\n";
	for (size_t i=0; i<g.getPatternGraph().getNodeNumber(); ++i) {
		out <<std::setw(6) <<i <<" ("<<g.getPatternGraph().getNodeLabel(i) <<")  --> ";
		for (sgm::Graph_Interface::OutEdge_iterator e = g.getPatternGraph().getOutEdgesBegin(i),
				eEnd = g.getPatternGraph().getOutEdgesEnd(i);
				e != eEnd; ++e)
		{
			out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
		}
		out <<" |\n";
	}
	if ( ! g.getConstraints().empty() ) {
		out <<" constraints : " <<g.getConstraints().size() <<"\n";
	}
	return out;
}
