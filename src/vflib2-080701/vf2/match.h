/*------------------------------------------------------------------
 * match.h
 * Header of match.cc
 * Declaration of the match function
 *
 * Author: P. Foggia
 * $Id: match.h,v 1.2 2011/11/04 12:35:30 mmann Exp $
 *-----------------------------------------------------------------*/



#ifndef MATCH_H
#define MATCH_H

#include "vf2/argraph.h"
#include "vf2/state.h"

namespace vf2 {

/*------------------------------------------------------------
 * Definition of the match_visitor type
 * a match visitor is a function that is invoked for
 * each match that has been found.
 * If the function returns false, then the next match is
 * searched; else the seach process terminates.
 -----------------------------------------------------------*/
typedef bool (*match_visitor)(int n, node_id c1[], node_id c2[], 
                              void *usr_data);

bool match(State *s0, int *pn, node_id c1[], node_id c2[]);

int match(State *s0, match_visitor vis, void *usr_data=NULL);

} // namespace vf2

#endif
