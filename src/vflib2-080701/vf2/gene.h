/*----------------------------------------------------
 * gene.h
 * Interface of gene.cc
 * Random generation of isomorphic ARGraphs
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: gene.h,v 1.2 2011/11/04 12:35:30 mmann Exp $
 ----------------------------------------------------*/

#ifndef GENE_H

#include "vf2/argraph.h"

namespace vf2 {

void Generate(int nodes, int edges, Graph **g1, Graph **g2, 
              bool connected=true);

} // namespace vf2

#endif
