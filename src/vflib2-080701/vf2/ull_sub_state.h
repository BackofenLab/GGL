/*------------------------------------------------------------------
 * ull_sub_state.h
 * Header of ull_sub_state.cc
 * Definition of the class UllSubState
 *
 * Author: P. Foggia
 * $Id: ull_sub_state.h,v 1.3 2011/11/04 12:35:29 mmann Exp $
 *-----------------------------------------------------------------*/


#ifndef ULL_SUB_STATE_H
#define ULL_SUB_STATE_H

#include "vf2/argraph.h"
#include "vf2/state.h"

namespace vf2 {

/*----------------------------------------------------------
 * class UllSubState
 * A representation of the current search state
 * of the Ullmann's algorithm for graph-subgraph isomorphism
 ---------------------------------------------------------*/
class UllSubState: public State
  { private:
	  node_id core_len;
      node_id *core_1;
      node_id *core_2;
      Graph *g1, *g2;
      node_id n1, n2;
      byte **M;   // Matrix encoding the compatibility of the nodes

      void refine();
    
    public:
      UllSubState(Graph *g1, Graph *g2);
      UllSubState(const UllSubState &state);
      ~UllSubState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1; };
      bool IsDead() { if (n1>n2) return true;
                      for(node_id i=core_len; i<n1; i++)
                        { for(node_id j=0; j<n2; j++)
                            if (M[i][j]!=0) goto next_row;
                          return true;
                      next_row: ;
                        }
                       return false;
                    };
                    node_id CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();
  };

} // namespace vf2

#endif
