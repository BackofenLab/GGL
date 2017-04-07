
#include "RRC_QuantumMechanics.hh"

#include <cassert>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <errno.h>
#include <cstdio>

// for server access
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	const std::string
	RRC_QuantumMechanics::
	defaultTCPhost = std::string("131.130.44.20");
	
	const int 
	RRC_QuantumMechanics::
	defaultTCPport = 9090;


////////////////////////////////////////////////////////////////////////////////

	
	RRC_QuantumMechanics
	::RRC_QuantumMechanics( const CalculationType calcType_
							, const std::string & tcpHost_
							, const int tcpPort_ )
	 :	calcType(calcType_)
	 	, tcpHost(tcpHost_)
	 	, tcpPort(tcpPort_)
	{
		// nothing else to do
	}

////////////////////////////////////////////////////////////////////////////////

	
	
	RRC_QuantumMechanics
	::~RRC_QuantumMechanics()
	{
		// nothing to do
	}

////////////////////////////////////////////////////////////////////////////////

	
	
	double
	RRC_QuantumMechanics
	::getRate( const ggl::chem::Reaction & reaction ) const
	{
		using namespace ggl;
		using namespace chem;

		double retVal = -1.0;

		// the string to create for the call
		// convert reaction into string: O.CCBr OC(C)Br OCC.Br MOPAC
		std::ostringstream call;

		/////////// COPY MOLECULES //////////////////////////////////////

		// add all molecules to call
		Reaction::Metabolite_Container::iterator mol = reaction.metabolites.begin();
		if (mol != reaction.metabolites.end()) {
			call << (*mol);
		}
		for (mol++; mol != reaction.metabolites.end(); ++mol) {
			call << "." << (*mol);
		}

		/////////// COPY TRANSITIONSTATE ////////////////////////////////

		call << " " << reaction.transState << " ";

		/////////// COPY PRODUCTS //////////////////////////////////////

		// add all products to call
		Reaction::Product_Container::iterator prod = reaction.products.begin();
		if (prod != reaction.products.end()) {
			call << (*prod);
		}
		for (prod++; prod != reaction.products.end(); ++prod) {
			call << "." << (*prod);
		}

		/////////// ADD CALCULATION TYPE ///////////////////////////////

		// set calculation type
		switch (calcType) {
		case QM_MOPAC:
			call << " MOPAC";
			break;
		case QM_JAGGUAR:
			call << " JAGGUAR";
			break;
		default:
			assert(false /*unhandled calculation type*/);
		}

		/////////// ACCESS SERVER //////////////////////////////////////


		// check if we can open a socket
		int sockfd = socket(AF_INET, SOCK_STREAM, 0);
		if (sockfd < 0) {
			std::cerr << "\n\tPROBLEM : Cannot open a socket : error code = "
					<< errno << "\n";
			return -1.0;
		}

		// set server details
		struct sockaddr_in serv_addr;
		memset(&serv_addr, 0, sizeof(serv_addr));
		serv_addr.sin_family = AF_INET;
		serv_addr.sin_port = htons(tcpPort);
		serv_addr.sin_addr.s_addr = inet_addr(tcpHost.c_str());

		// try connection
		if (connect(sockfd, (struct sockaddr*) &serv_addr, sizeof(serv_addr)) < 0) {
			std::cerr << "\n\tPROBLEM : Cannot connect to server : error code = "
					<< errno << "\n";
			return -1.0;
		}

		// write the call to server
		if (write(sockfd, call.str().c_str(), call.str().size()) < 0) {
			std::cerr << "\n\tPROBLEM : Cannot write to socket : error code = "
					<< errno << "\n";
			return -1.0;
		}

		// read the answer from server
		char rateAnswerBuffer[30];
		memset(rateAnswerBuffer, 0, sizeof(rateAnswerBuffer));
		if (read(sockfd, rateAnswerBuffer, sizeof(rateAnswerBuffer)) < 0) {
			std::cerr << "\n\tPROBLEM : Cannot read socket : error code = "
					<< errno << "\n";
			return -1.0;
		}

		// parse return value from answer buffer
		sscanf(rateAnswerBuffer, "%lf", &retVal);

		// return calculated rate
		return retVal;
	}

////////////////////////////////////////////////////////////////////////////////

	
	
	bool
	RRC_QuantumMechanics
	::needTransitionState( void ) const
	{
		  // transition state needed for this rates .. 
		return true;
	}

////////////////////////////////////////////////////////////////////////////////

	int
	RRC_QuantumMechanics::
	testServerAvailability( void ) const
	{
	    // the string to create for the call
	    // convert reaction into string: O.CCBr OC(C)Br OCC.Br MOPAC
		std::ostringstream call;
		
		  // dummy call
		call <<"O.CCBr OC(C)Br OCC.Br";

	    // set calculation type
		switch (calcType) {
		case QM_MOPAC :	call <<" MOPAC"; break;
		case QM_JAGGUAR :	call <<" JAGGUAR"; break;
		default :
			assert( false /*unhandled calculation type*/ );
		}

		
		/////////// ACCESS SERVER //////////////////////////////////////


		    // check if we can open a socket
		  int sockfd = socket(AF_INET, SOCK_STREAM, 0);
		  if (sockfd < 0) {
			  return errno;
		  }

		    // set server details
		  struct sockaddr_in serv_addr;
		  memset(&serv_addr, 0, sizeof(serv_addr));
		  serv_addr.sin_family = AF_INET;
		  serv_addr.sin_port = htons(tcpPort);
		  serv_addr.sin_addr.s_addr = inet_addr(tcpHost.c_str());

		    // try connection
		  if (connect(sockfd, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
			  return errno;
		  }

		    // write the call to server
		  if (write(sockfd, call.str().c_str(), call.str().size()) < 0) {
			  return errno;
		  }

		    // read the answer from server
		  char rateAnswerBuffer[30];
		  memset(rateAnswerBuffer, 0, sizeof(rateAnswerBuffer));
		  if (read(sockfd, rateAnswerBuffer, sizeof(rateAnswerBuffer)) < 0) {
			  return errno;
		  }

		    // everything ok ..
		  return 0;
	}

	
////////////////////////////////////////////////////////////////////////////////

	std::string
	RRC_QuantumMechanics::
	getTCPhost( void ) const {
		return tcpHost;
	}
	
	
////////////////////////////////////////////////////////////////////////////////

	int
	RRC_QuantumMechanics::
	getTCPport( void ) const {
		return tcpPort;
	}

	
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

