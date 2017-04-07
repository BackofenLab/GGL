/*-----------------------------------------------------
 * error.cc
 * Implementation of error handling functions
 *
 * Author: P. Foggia
 * $Id: error.cc,v 1.4 2011/11/04 12:35:29 mmann Exp $
 ----------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "vf2/error.h"

namespace vf2 {

/*------------------------------------------
 * void error(msg, ...)
 * Prints an error message and exits 
 * the program.
 * The syntax is the same of printf, 
 * except that a trailing \n is automatically
 * appended.
 -----------------------------------------*/
void error(const char *msg, ...)
  { va_list ap;
    va_start(ap, msg);
    fprintf(stderr, "ERROR: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
  }

} // namespace vf2
