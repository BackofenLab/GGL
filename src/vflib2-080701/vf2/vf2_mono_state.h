/*------------------------------------------------------------
 * vf2_mono_state.h
 * Interface of vf2_mono_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs.
 * See: argraph.h state.h
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/




#ifndef VF2_MONO_STATE_H
#define VF2_MONO_STATE_H

#include "vf2/argraph.h"
#include "vf2/state.h"

namespace vf2 {

/*----------------------------------------------------------
 * class VF2MonoState
 * A representation of the SSR current state
 * See vf2_sub_state.cc for more details.
 ---------------------------------------------------------*/
class VF2MonoState: public State
  { typedef ARGraph_impl Graph;

    private:
    	node_id core_len, orig_core_len;
    	node_id added_node1;
      node_id t1both_len, t2both_len, t1in_len, t1out_len, 
          t2in_len, t2out_len; // Core nodes are also counted by these...
      node_id *core_1;
      node_id *core_2;
      node_id *in_1;
      node_id *in_2;
      node_id *out_1;
      node_id *out_2;

      node_id *order;

      Graph *g1, *g2;
      node_id n1, n2;

	  long *share_count;
    
    public:
      VF2MonoState(Graph *g1, Graph *g2, bool sortNodes=false);
      VF2MonoState(const VF2MonoState &state);
      ~VF2MonoState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1 ; };
      bool IsDead() { return n1>n2  || 
                         t1both_len>t2both_len ||
                         t1out_len>t2out_len ||
                         t1in_len>t2in_len;
                    };
        node_id CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();

      virtual void BackTrack();
  };

} // namespace vf2

#endif

