/*----------------------------------------------------
 * gene_mesh.h
 * Interface of gene_mesh.cc
 * Random generation of isomorphic ARGraphs which are
 * almost meshes.
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: gene_mesh.h,v 1.2 2011/11/04 12:35:31 mmann Exp $
 ----------------------------------------------------*/

#ifndef GENE_MESH_H

#include "vf2/argraph.h"

namespace vf2 {

void GenerateMesh(int nodes, int extra_edges, Graph **g1, Graph **g2,
                  int sub_nodes=-1);

} // namespace vf2

#endif
