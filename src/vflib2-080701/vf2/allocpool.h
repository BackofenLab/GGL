/*-----------------------------------------------------
 * allocpool.h
 * Definition of an allocation pool class.
 *
 * Author: P. Foggia
 * $Id: allocpool.h,v 1.2 2011/11/04 12:35:29 mmann Exp $
 ----------------------------------------------------*/

/*----------------------------------------------------------------------
 *   CLASS DESCRIPTION
 * An allocation pool is an efficient way to allocate on heap
 * a set of objects of the same type that will be deleted
 * at the same time. Instead of allocating each object separately,
 * thus incurring in both space and time overhead, "chunks" of 
 * objects are allocated. The allocation pool mantains a list of
 * the chunks, that will be deleted by the pool distructor.
 *
 * Two classes are provided: an Allocator, which is an
 * interface independent of actual allocation strategy, and 
 * AllocationPool which is the real class. 
 ---------------------------------------------------------------------*/

#ifndef ALLOCPOOL_H
#define ALLOCPOOL_H


#include <stddef.h>
#include "error.h"

namespace vf2 {

/*-----------------------------------------------------------
 * Declaration of class Allocator,
 * which represents an abstract way of allocating objects
 *---------------------------------------------------------*/

template <class T>
class Allocator
  { protected:
      virtual T *alloc() = 0;
    public:
      virtual ~Allocator() {}
      virtual T *Allocate() 
          { T*p=alloc();
            if (p==NULL)
              error("Out of memory");
            return p;
          }
  };


/*-----------------------------------------------------------
 * Declaration of class NewAllocator,
 * which simply use new to allocate objects
 ----------------------------------------------------------*/
template <class T>
class NewAllocator: public Allocator<T>
  { protected:
      virtual T* alloc() { return new T; }
  };

/*-----------------------------------------------------------
 * Declaration of class NullAllocator, which always
 * return a null pointer
 ----------------------------------------------------------*/
template <class T>
class NullAllocator: public Allocator<T>
  { public:
      virtual T *Allocate() { return NULL; }
	protected:
	  virtual T *alloc() { return NULL; }
  };


/*-----------------------------------------------------------
 * Declaration of class AllocationPool
 *---------------------------------------------------------*/

template <class T, int CHUNK_SIZE>
class AllocationPool: public Allocator<T>
  { private:
      struct chunk
        { chunk *link;
          T content[CHUNK_SIZE];
        };
      chunk *chunkList;
      int rest;
      void grow();
	protected:
      virtual T *alloc();
    public:
      AllocationPool();
      ~AllocationPool();
  };



/*----------------------------------------------------
 * inline methods of class AllocationPool
 ---------------------------------------------------*/
template <class T, int CHUNK_SIZE>
AllocationPool<T, CHUNK_SIZE>::AllocationPool() 
  { chunkList=0; 
    rest=0; 
  }

template <class T, int CHUNK_SIZE>
AllocationPool<T, CHUNK_SIZE>::~AllocationPool() 
  { chunk *p=chunkList;
    while (p!=0)
      { chunkList=p->link;
        delete p;
        p=chunkList;
      }
  }

template <class T, int CHUNK_SIZE>
T * AllocationPool<T, CHUNK_SIZE>::alloc() 
  { if (rest==0)
      grow();
    if (rest>0)
      return chunkList->content + (--rest);
    else
      return 0;
  }

template <class T, int CHUNK_SIZE>
void AllocationPool<T, CHUNK_SIZE>::grow() 
  { chunk *p=new chunk;
    if (p!=0)
      { p->link=chunkList;
        chunkList=p;
        rest=CHUNK_SIZE;
      }
  }

}

#endif
