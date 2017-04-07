/*------------------------------------------------------------
 * state.h
 * Definition of an abstract class representing a state of the 
 * matching process between two ARGs.
 * See: argraph.h 
 *
 * Author: P. Foggia
 * $Id: state.h,v 1.3 2011/11/04 12:35:29 mmann Exp $
 *-----------------------------------------------------------------*/



#ifndef STATE_H
#define STATE_H

#include "vf2/argraph.h"

namespace vf2 {

/*----------------------------------------------------------
 * class State
 * An abstract representation of the SSR current state.
 * NOTE: Respect to pre-2.0 version of the library, class
 *   State assumes explicitly a depth-first search. The
 *   BackTrack method has been added to allow a state 
 *   to clean up things before reverting to its parent. This
 *   can be used, for instance, for sharing resources between
 *   the parent and the child. The BackTrack implementation
 *   can safely assume that at most one AddPair has been
 *   performed on the state.
 ---------------------------------------------------------*/
class State
  { 

    public:
      virtual ~State() {} 
      virtual Graph *GetGraph1()=0;
      virtual Graph *GetGraph2()=0;
      virtual bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE)=0;
      virtual bool IsFeasiblePair(node_id n1, node_id n2)=0;
      virtual void AddPair(node_id n1, node_id n2)=0;
      virtual bool IsGoal() =0;
      virtual bool IsDead() =0;
      virtual node_id CoreLen() =0;
      virtual void GetCoreSet(node_id c1[], node_id c2[]) =0;
      virtual State *Clone() =0;  // Changed clone to Clone for uniformity
     
      virtual void BackTrack() { };
  };

} // namespace vf2

#endif

