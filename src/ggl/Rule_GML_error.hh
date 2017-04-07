/*
 * Rule_GML_error.hh
 *
 *  Created on: 29.09.2010
 *      Author: mmann
 */

#ifndef GGL_RULE_GML_ERROR_HH_
#define GGL_RULE_GML_ERROR_HH_

#include <stdexcept>
#include <string>


namespace ggl {

	/*! @brief Rule parsing error
	 *
	 * Specialized error class to report GML rule parsing errors like multiple
	 * use of indices etc.
	 *
	 * @author Martin Mann (c) 2010 http://www.bioinf.uni-freiburg.de/~mmann/
	 *
	 */
    class Rule_GML_error : public std::logic_error
    {
    public:
      explicit Rule_GML_error (const std::string& what_arg)
		  : std::logic_error( what_arg )
      {}
	};


} // namespace gml

#endif /* GGL_RULE_GML_ERROR_HH_ */
