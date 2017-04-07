/*------------------------------------------------------------
 * vf_mono_state.h
 * Interface of vf_mono_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs, for the monomorphism problem
 * See: argraph.h state.h
 *
 * Author: P. Foggia
 * $Id: vf_mono_state.h,v 1.3 2011/11/04 12:35:30 mmann Exp $
 *-----------------------------------------------------------------*/



#ifndef VF_MONO_STATE_H
#define VF_MONO_STATE_H

#include "vf2/argraph.h"
#include "vf2/state.h"


namespace vf2 {

/*----------------------------------------------------------
 * class VFMonoState
 * A representation of the SSR current state
 * for a graph-subgraph isomorphism.
 ---------------------------------------------------------*/
class VFMonoState: public State
  { typedef ARGraph_impl Graph;

    private:
      // Flags used to encode the state
      enum { ST_CORE=0x01,
             ST_TERM_IN=0x02,
             ST_TERM_OUT=0x04
            };

    node_id core_len;
    node_id t1in_len, t1out_len, t2in_len, t2out_len;
      node_id *core_1;
      node_id *core_2;
      byte *node_flags_1;
      byte *node_flags_2;
      Graph *g1, *g2;
      node_id n1, n2;
    
    public:
      VFMonoState(Graph *g1, Graph *g2);
      VFMonoState(const VFMonoState &state);
      ~VFMonoState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1; };
      bool IsDead() { return n1>n2 || 
                      (t1out_len>t2out_len) ||
                      (t1in_len>t2in_len); 
                    };
    node_id CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();
  };

} // namespace vf2

#endif

