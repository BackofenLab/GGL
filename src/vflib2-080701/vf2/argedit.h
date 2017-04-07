/*------------------------------------------------------------------
 * argedit.h
 * Definition of a simple ARG loader which allows graph edit
 * operations
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: argedit.h,v 1.3 2011/11/04 12:35:29 mmann Exp $
 *-----------------------------------------------------------------*/


#ifndef ARGEDIT_H
#define ARGEDIT_H

#include "vf2/argraph.h"

namespace vf2 {

/*---------------------------------------------------------
 * Class ARGEdit
 * A simple ARGLoader providing graph edit operations.
 * Note: the ARGEdit does not make provisions for the
 *       deallocation of the attributes, which must be
 *       dealt with by the programmer.
 -------------------------------------------------------*/
class ARGEdit: public ARGLoader
  { 
    public:

      ARGEdit();
      ARGEdit(ARGraph_impl &g);
      ARGEdit(ARGLoader &g);
      ~ARGEdit();

      /* Redefined ARGLoader methods */
      virtual node_id NodeCount();
      virtual void *GetNodeAttr(node_id node);
      virtual SIZETYPE OutEdgeCount(node_id node);
      virtual node_id GetOutEdge(node_id node, SIZETYPE i, void **pattr);

      /* Graph edit operations */
      node_id InsertNode(void *attr);
      void InsertEdge(node_id n1, node_id n2, void *attr);
      void DeleteNode(node_id n);
      void DeleteEdge(node_id n1, node_id n2);

    protected:
    	node_id count;

      struct eNode
        { node_id from;
          node_id to;
          SIZETYPE pos;
          void *attr;
          eNode *next;
        };

      struct nNode 
        { node_id id;
          SIZETYPE count;
          void *attr;
          nNode *next;
          eNode *edges;
        };

      nNode *nodes;
      nNode *lastNode;
      eNode *lastEdge;

      virtual void destroyNodeAttr(void *) {};
      virtual void destroyEdgeAttr(void *) {};
  };
  
} // namespace vf2

#endif
