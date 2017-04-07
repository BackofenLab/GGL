#ifndef GGL_CHEM_REACTION_HH_
#define GGL_CHEM_REACTION_HH_


#include <string>
#include <set>
#include <ostream>

namespace ggl {
 namespace chem {

	  /*! @brief Reaction description
	   *
	   * A data container describing information of a chemical reaction.
	   * 
	   * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	   * @author Alexander Ullrich
	   * 
	   */
	class Reaction  {
		
	public:
		
		  //! container that stores all metabolites/educts of the reaction
		typedef std::multiset< std::string > Metabolite_Container;

		  //! container that stores all products of the reaction
		typedef std::multiset< std::string > Product_Container;
		
	public:

		  /*! the identifier of the ggl::chem::ChemRule rule that was applied
		   * in this reaction.
		   */
		std::string rule_id;
		
		  //! the SMILES of the molecules involved in this reaction
		Metabolite_Container metabolites;
		
		  //! the SMILES of the produced molecules of this reaction
		Product_Container products;
		
		  //! the rate of the reaction
		double rate;
		
		  //! The SMILES of the transition state along the reaction. It is used
		  //! to calculate appropriate reaction rates.
		std::string transState;
		
	public:
		
		  //! Default construction with empty content
		Reaction();
		
		  /*! Construction
		   * @param rule_id the identifier of the rule that was applied
		   * @param metabolites the SMILES of the molecules involved
		   * @param products the SMILES of the produced molecules
		   * @param rate the rate of the reaction
		   * @param transState the SMILES of the transition state
		   */
		Reaction(	const std::string& rule_id
					, const Metabolite_Container& metabolites
					, const Product_Container& products
					, const double rate
					, const std::string& transState );

		//! Destruction
		~Reaction();
		
		  //! resets the content to values of default construction;
		void
		clear( void );
		
		  /*! Comparison operator to enable the storage of Reaction objects in
		   * ordered containers like std::set that require a strict less
		   * ordering.
		   * 
		   * @param r2 the Reaction to compare to
		   * @return true if [
		   *      (rule_id < r2.rule_id) or equal and
		   *      (less products) or equal and
		   *      (less metabolites) or equal and
		   *      (equal product number and one SMILES smaller) or equal and
		   *      (equal metabolites number and one SMILES smaller) or equal and
		   *      (both rates are set) and (rate < r2.rate) ],
		   *      false otherwise
		   */
		bool
		operator < ( const Reaction& r2 ) const ;
		
		
	};

 } // namespace chem
} // namespace ggl
	
	 /*! outstream operator for Reaction class that prints the reaction to stream
	  * @param out the stream to write to
	  * @param r the Reaction to write
	  * @return the modified stream
	  */
	std::ostream&
	operator << ( std::ostream& out, const ggl::chem::Reaction& r );

	
#endif /*REACTION_HH_*/

		
