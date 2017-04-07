#ifndef RRC_QUANTUMMECHANICS_HH_
#define RRC_QUANTUMMECHANICS_HH_


#include "ggl/chem/ReactionRateCalculation.hh"

#include <string>


	 /*! Computes a reaction rate based on the transition state of an reaction.
	  * It utilizes an according server written by Sebastian Sauer.
	  * 
	  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mann/
	  * @author Daniel Hoegerl
	  */
	class RRC_QuantumMechanics : public ggl::chem::ReactionRateCalculation
	{
	public:
		
		  //! default server URL
		static const std::string defaultTCPhost;
		
		  //! default server port
		static const int defaultTCPport;
		
		  //! Available options for reaction rate calculation
		enum CalculationType { QM_MOPAC, QM_JAGGUAR };
		
	protected:
		
		  //! calculation type
		CalculationType calcType;
		
		  //! server URL or IP to use for access
		std::string tcpHost;
		
		  //! server PORT to use for access
		int tcpPort;
		
		
	public:
		
		
		 /*! Construction of the rate calculation object that sets the server
		  * details.
		  * 
		  * @param tcpHost the URL of the server
		  * @param tcpPort the port to reach the server
		  */
		RRC_QuantumMechanics(	const CalculationType = QM_MOPAC
								, const std::string & tcpHost = defaultTCPhost
								, const int tcpPort = defaultTCPport );
		
		
		virtual
		~RRC_QuantumMechanics();
		
		  /*! Calculates the reaction rate for a given Reaction via a call to
		   * the reaction rate calculation server.
		   * 
		   * @param reaction the Reaction object to calculate the rate for
		   * @return the according reaction rate or -1.0 if an error occured. a
		   *         corresponding message is written to std::cerr.
		   */
		virtual
		double
		getRate( const ggl::chem::Reaction & reaction ) const;
		
		  /*! Announces that the reaction rate calculation needs the explicit 
		   * transition state within the reaction information.
		   * 
		   * @return true 
		   */
		virtual
		bool
		needTransitionState( void ) const;

		
		  /*! Tries to access the server and submits a dummy request.
		   * 
		   * @return 0 if the server is working properly; an error code != 0
		   *         otherwise.
		   */
		int
		testServerAvailability( void ) const;
		
		  /*! Access to the set TCP host.
		   * @return the currently set TCP host URL.
		   */
		std::string
		getTCPhost( void ) const;
		
		  /*! Access to the set TCP port.
		   * @return the currently set TCP port.
		   */
		int
		getTCPport( void ) const;
		
	};


#endif /*RRC_QUANTUMMECHANICS_HH_*/
