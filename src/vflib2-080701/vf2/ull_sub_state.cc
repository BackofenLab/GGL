/*------------------------------------------------------------
 * ull_sub_state.cc
 * Implementation of the class UllSubState
 *
 * Author: P. Foggia
 * $Id: ull_sub_state.cc,v 1.4 2011/11/19 23:05:53 mmann Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * NOTES: 
 *   g1 is the graph matched with a subgraph of g2.
 *   The attribute compatibility check (methods CompatibleNode
 *   and CompatibleEdge of ARGraph) is always performed
 *   applying the method to g1, and passing the attribute of
 *   g1 as first argument, and the attribute of g2 as second
 *   argument. This may be important if the compatibility
 *   criterion is not symmetric.
 -----------------------------------------------------------------*/


#include "vf2/ull_sub_state.h"

#include "vf2/error.h"

namespace vf2 {

/*----------------------------------------------------------
 * Methods of the class UllSubState
 ---------------------------------------------------------*/

/*----------------------------------------------------------
 * UllSubState::UllSubState(g1, g2)
 * Constructor. Makes an empty state.
 ---------------------------------------------------------*/
UllSubState::UllSubState(Graph *ag1, Graph *ag2)
  { g1=ag1;
    g2=ag2;
    n1=g1->NodeCount();
    n2=g2->NodeCount();

    core_len=0;

    core_1=new node_id[n1];
    core_2=new node_id[n2];
    M=new byte *[n1];
    if (!core_1 || !core_2 || !M)
      error("Out of memory");

    node_id i,j;

    for (i=0; i<n1; ++i )
      { M[i]=new byte[n2];
        if (!M[i])
          error("Out of memory");
      }

    for(i=0; i<n1; ++i )
      { 
        core_1[i]=NULL_NODE;
      }
    for(i=0; i<n2; ++i )
      {
        core_2[i]=NULL_NODE;
      }
    for(i=0; i<n1; ++i )
      for(j=0; j<n2; j++)
        M[i][j]=(g1->InEdgeCount(i) <= g2->InEdgeCount(j) &&
                 g1->OutEdgeCount(i) <= g2->OutEdgeCount(j)) &&
                 g1->CompatibleNode(g1->GetNodeAttr(i), g2->GetNodeAttr(j)) ?
                 1: 0;

  }


/*----------------------------------------------------------
 * UllSubState::UllSubState(state)
 * Copy constructor. 
 ---------------------------------------------------------*/
UllSubState::UllSubState(const UllSubState &state)
  { g1=state.g1;
    g2=state.g2;
    n1=state.n1;
    n2=state.n2;

    core_len=state.core_len;

    core_1=new node_id[n1];
    core_2=new node_id[n2];
    M=new byte *[n1];
    if (!core_1 || !core_2 || !M)
      error("Out of memory");

    node_id i,j;

    for (i=0; i<core_len; ++i )
      M[i]=NULL;
 
    for (i=core_len; i<n1; ++i )
      { M[i]=new byte[n2];
        if (!M[i])
          error("Out of memory");
      }

    for(i=0; i<n1; ++i )
      core_1[i]=state.core_1[i];
    for(i=0; i<n2; ++i )
        core_2[i]=state.core_2[i];
    for(i=core_len; i<n1; ++i )
      for(j=0; j<n2; j++)
        M[i][j]=state.M[i][j];
  }


/*---------------------------------------------------------------
 * UllSubState::~UllSubState()
 * Destructor.
 --------------------------------------------------------------*/
UllSubState::~UllSubState() 
  { delete [] core_1;
    delete [] core_2;
    node_id i;
    for(i=0; i<n1; ++i )
      if (M[i])
        delete [] M[i];
    delete [] M;
  }


/*--------------------------------------------------------------------------
 * bool UllSubState::NextPair(pn1, pn2, prev_n1, prev_n2)
 * Puts in *pn1, *pn2 the next pair of nodes to be tried.
 * prev_n1 and prev_n2 must be the last nodes, or NULL_NODE (default)
 * to start from the first pair.
 * Returns false if no more pairs are available.
 -------------------------------------------------------------------------*/
bool UllSubState::NextPair(node_id *pn1, node_id *pn2,
              node_id prev_n1, node_id prev_n2)
  { if (prev_n1==NULL_NODE)
      { prev_n1=core_len;
        prev_n2=0;
      }
    else if (prev_n2==NULL_NODE)
      prev_n2=0;
    else
      prev_n2++;

    if (prev_n2>=n2)
      { prev_n1++;
        prev_n2=0;
      }

    if (prev_n1!=core_len)
      return false;
    while (prev_n2<n2 && M[prev_n1][prev_n2]==0)
      prev_n2++;
    if (prev_n2<n2)
      { *pn1=prev_n1;
        *pn2=prev_n2;
        return true;
      }
    else
      return false; 
  }



/*---------------------------------------------------------------
 * bool UllSubState::IsFeasiblePair(node1, node2)
 * Returns true if (node1, node2) can be added to the state
 --------------------------------------------------------------*/
bool UllSubState::IsFeasiblePair(node_id node1, node_id node2)
  { assert(node1<n1);
    assert(node2<n2);

    return M[node1][node2]!=0;
  }



/*--------------------------------------------------------------
 * void UllSubState::AddPair(node1, node2)
 * Adds a pair to the Core set of the state.
 * Precondition: the pair must be feasible
 -------------------------------------------------------------*/
void UllSubState::AddPair(node_id node1, node_id node2)
  { assert(node1<n1);
    assert(node2<n2);
    assert(core_len<n1);
    assert(core_len<n2);

    core_1[node1]=node2;
    core_2[node2]=node1;

    core_len++;

    node_id k;

    for(k=core_len; k<n1; k++)
      M[k][node2]=0;

    refine(); 
  }



/*--------------------------------------------------------------
 * void UllSubState::GetCoreSet(c1, c2)
 * Reads the core set of the state into the arrays c1 and c2.
 * The i-th pair of the mapping is (c1[i], c2[i])
 --------------------------------------------------------------*/
void UllSubState::GetCoreSet(node_id c1[], node_id c2[])
  { node_id i,j;
    for (i=0,j=0; i<n1; ++i )
      if (core_1[i] != NULL_NODE)
        { c1[j]=i;
          c2[j]=core_1[i];
          j++;
        }
  }



/*------------------------------------------------------------
 * void UllSubState::refine()                             PRIVATE
 * Removes from the matrix M all the pairs which are not
 * compatible with the isomorphism condition
 -----------------------------------------------------------*/
void UllSubState::refine()
  {
	node_id i, j, k, l;

    for(i=core_len; i<n1; ++i )
      for(j=0; j<n2; j++)
        if (M[i][j])
          { bool edge_ik, edge_ki, edge_jl, edge_lj;
		    // The following (commented-out) for wasn't necessary
            //for(k=0; k<core_len; k++)
			for(k=core_len-1; k<core_len; k++)
              { l=core_1[k];
                assert(l!=NULL_NODE);
                edge_ik=g1->HasEdge(i,k);
                edge_ki=g1->HasEdge(k,i);
                edge_jl=g2->HasEdge(j,l);
                edge_lj=g2->HasEdge(l,j);
                if (edge_ik!=edge_jl || edge_ki!=edge_lj)
                  { M[i][j]=0;
                    break;
                  }
                else if (edge_ik  &&
                         !g1->CompatibleEdge(g1->GetEdgeAttr(i,k),
                                             g2->GetEdgeAttr(j,l)))
                  { M[i][j]=0;
                    break;
                  }
                else if (edge_ki  &&
                         !g1->CompatibleEdge(g1->GetEdgeAttr(k,i),
                                             g2->GetEdgeAttr(l,j)))
                  { M[i][j]=0;
                    break;
                  }
              }
          }
  }


/*-----------------------------------------
 * Clones a state, allocating with new
 ----------------------------------------*/
State* UllSubState::Clone()
  { return new UllSubState(*this);
  }

} // namespace vf2

