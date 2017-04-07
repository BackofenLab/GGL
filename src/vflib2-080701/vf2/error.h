/*----------------------------------------------------
 * error.h
 * Header of error.cc
 * Error handling functions
 *
 * Author: P. Foggia
 * $Id: error.h,v 1.4 2011/11/04 12:35:30 mmann Exp $
 *--------------------------------------------------*/

#ifndef ERROR_H

#include <stddef.h>

namespace vf2 {

void error(const char *msg, ...);

#define FAIL(reason)    error("%s in %s:%d", (reason), __FILE__, \
                                                          (int)__LINE__)

#define OUT_OF_MEMORY()  FAIL("Out of memory")
#define CANT_HAPPEN()    FAIL("Can't happen")



} // namespace vf2

#endif
