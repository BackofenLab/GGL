
#include "ggl/chem/Reaction.hh"
#include <limits>


namespace ggl {
 namespace chem {


	  // Default construction with empty content
	Reaction::
	Reaction()
	 :	rule_id("")
	 	, metabolites()
	 	, products()
	 	, rate(std::numeric_limits<double>::signaling_NaN())
	 	, transState()
	{
		clear(); // to ensure the same values when calling clear later
	}
		
	  /* Construction
	   * @param rule_id_ the identifier of the rule that was applied
	   * @param metabolites_ the SMILES of the molecules involved
	   * @param products_ the SMILES of the produced molecules
	   * @param rate_ the rate of the reaction
	   */
	Reaction::
	Reaction(	const std::string& rule_id_
				, const Metabolite_Container& metabolites_ 
				, const Product_Container& products_ 
				, const double rate_
				, const std::string& transState_)
	 :	rule_id( rule_id_ )
	 	, metabolites( metabolites_ )
	 	, products( products_ )
	 	, rate( rate_ )
	 	, transState( transState_ )
	{}

	// Destruction
	Reaction::
	~Reaction()
	{}
		
	  // resets the content to values of default construction;
	void
	Reaction::
	clear( void )
	{
		rule_id = std::string("");
	 	metabolites.clear();
	 	products.clear();
	 	rate = std::numeric_limits<double>::signaling_NaN();
	 	transState = std::string("");
	}
		
	  /* Comparison operator to enable the storage of Reaction objects in
	   * ordered containers like std::set that require a strict less
	   * ordering.
	   * 
	   * @param r2 the Reaction to compare to
	   * @return true if 
	   *      (rule_id < r2.rule_id) or equal and
	   *      (less products) or equal and
	   *      (less metabolites) or equal and
	   *      (equal product number and one SMILES smaller) or equal and
	   *      (equal metabolites number and one SMILES smaller) or equal and
	   *      (rate < r2.rate),
	   *      false otherwise
	   */
	bool
	Reaction::
	operator < ( const Reaction& r2 ) const 
	{
		 // compare rule_id
		int cmp = rule_id.compare( r2.rule_id );
		if (cmp < 0) 
			return true;
		if (cmp > 0)
			return false;

		 // compare metabolite number
		if (this->metabolites.size() < r2.metabolites.size())
			return true;
		if (this->metabolites.size() > r2.metabolites.size())
			return false;

		 // compare product number
		if (this->products.size() < r2.products.size())
			return true;
		if (this->products.size() > r2.products.size())
			return false;
		
		 // compare metabolites directly
		Metabolite_Container::const_iterator m1 = this->metabolites.begin();
		Metabolite_Container::const_iterator m2 = r2.metabolites.begin();
		while (m1 != this->metabolites.end()) {
			cmp = m1->compare(*m2);
			if (cmp < 0) return true;
			if (cmp > 0) return false;
		    m1++; m2++;
		}
		
		 // compare transition state string
		cmp = this->transState.compare(r2.transState);
		if (cmp < 0) return true;
		if (cmp > 0) return false;

		 // compare products directly
		Product_Container::const_iterator p1 = this->products.begin();
		Product_Container::const_iterator p2 = r2.products.begin();
		while (p1 != this->products.end()) {
			cmp = p1->compare(*p2);
			if (cmp < 0) return true;
			if (cmp > 0) return false;
		    p1++; p2++;
		}
		
		
		 // check if both rates are not NaN and compare rates
		return	(this->rate == this->rate)	// false if NaN
				&& (r2.rate == r2.rate)		// false if NaN
				&& (this->rate < r2.rate);
	}
		
		

 } // namespace chem
} // namespace ggl

	
	 /* outstream operator for Reaction class that prints the reaction to stream
	  * @param out the stream to write to
	  * @param r the Reaction to write
	  * @return the modified stream
	  */
	std::ostream&
	operator << ( std::ostream& out, const ggl::chem::Reaction& r )
	{
		ggl::chem::Reaction::Metabolite_Container::const_iterator m = r.metabolites.begin();
		ggl::chem::Reaction::Product_Container::const_iterator p = r.products.begin();
		
		 // write rule_id and first of metabolites
		out <<r.rule_id 
			<<" : " 
			<<(r.metabolites.size() > 0 ? (*m) : "") ;
		 // write remaining metabolites
		if ( r.metabolites.size() > 0) {
			for ( m++; m != r.metabolites.end(); ++m)
				out <<"." <<(*m);
		}
		
		 // write reaction separator and transition state if present
		out <<">" <<r.transState <<">";
		
		 // write first of products
		out	<<(r.products.size() > 0 ? (*p) : "") ;
		 // write remaining products
		if ( r.products.size() > 0 ) {
			for ( p++; p != r.products.end(); ++p)
				out <<"." <<(*p);
		}
		 // write reaction rate
		out <<" : " 
			<<r.rate;
		
		return out;
	}

