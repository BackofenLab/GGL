#ifndef GGL_CHEM_REACTIONRATECALCULATION_HH_
#define GGL_CHEM_REACTIONRATECALCULATION_HH_

#include "ggl/chem/Reaction.hh"

#include <ggl/Rule.hh>


namespace ggl {
 namespace chem {

	////////////////////////////////////////////////////////////////////////////

	 /*! @brief Interface reaction rate calculation
	  *
	  * Abstract super class for reaction rate calculations.
	  * 
	  * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	  * 
	  */
	class ReactionRateCalculation
	{
	public:
		
		  //! default construction
		ReactionRateCalculation() {}
		
		  //! default destruction
		virtual ~ReactionRateCalculation() {}
		
		  /*! Calculates the reaction rate for a given Reaction.
		   * 
		   * @param reaction the Reaction object to calculate the rate for
		   * @return the according reaction rate
		   */
		virtual
		double
		getRate( const Reaction & reaction ) const = 0;
		
		  /*! Announces whether the reaction rate calculation needs the explicit 
		   * transition state present within the reaction information or not.
		   * 
		   * @return true iff the reaction calculation needs the presence of an
		   *         explicit transition state, false otherwise (no transition
		   *         state needed)
		   */
		virtual
		bool
		needTransitionState( void ) const = 0;
	};

	////////////////////////////////////////////////////////////////////////////


	 /*! A reaction rate calculator that does not calculate a rate but returns
	  * a fixed reaction rate independently from the given reaction.
	  * 
	  * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	  * 
	  */
	class RRC_Fixed : public ReactionRateCalculation
	{
	protected:
		
		  //! the fixed rate to return for all reactions
		const double fixRate;
		
	public:
		
		  /*! Constructs a new reaction rate calculator that yields always the
		   * same fixed reaction rate given.
		   * 
		   * @param fixRate the fixed rate to be assigned to all reactions
		   */
		RRC_Fixed(const double fixRate );
		
		  //! default destruction
		virtual ~RRC_Fixed();
		
		  /*! Returns a fixed transition rate for all reactions.
		   * 
		   * @param reaction the Reaction object to calculate the rate for; but
		   *          doesnt matter because fixed rate returned!
		   * @return the fixed reaction rate
		   */
		virtual
		double
		getRate( const Reaction & reaction ) const;
		
		  /*! No transition state needed, such that this returns always false!
		   * 
		   * @return false because no transition state needed for a fixed 
		   *          reaction rate
		   */
		virtual
		bool
		needTransitionState( void ) const;
	};


	////////////////////////////////////////////////////////////////////////////


	 /*! A reaction rate calculator that does not calculate a rate (returns 
	  * always "nan") but enforces, that the reactions transition state is
	  * calculated.
	  * 
	  * @author Martin Mann (c) 2009 http://www.bioinf.uni-freiburg.de/~mmann/
	  * 
	  */
	class RRC_TState : public ReactionRateCalculation
	{
	public:
		
		  /*! Constructs a new reaction rate calculator that yields always the
		   * same reaction rate ("nan") but enforces the creation of the 
		   * reactions transition state.
		   */
		RRC_TState();
		
		  //! default destruction
		virtual ~RRC_TState();
		
		  /*! Returns a "nan" transition rate for all reactions.
		   * 
		   * @param reaction the Reaction object to calculate the rate for; but
		   *          doesnt matter because "nan" returned!
		   * @return "nan"
		   */
		virtual
		double
		getRate( const Reaction & reaction ) const;
		
		  /*! YES, a transition state is needed!
		   * 
		   * @return true, because that's the class for! ;)
		   */
		virtual
		bool
		needTransitionState( void ) const;
	};


	////////////////////////////////////////////////////////////////////////////

 } // namespace chem
} // namespace ggl
	
#endif /*REACTIONRATECALCULATION_HH_*/

