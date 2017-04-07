/*--------------------------------------------------------
 * xsubgraph.h
 * Interface of xsubgraph.cc
 * Random extraction of a (possibly) connected subgraph
 * See: argraph.h
 *
 * Author: P. Foggia
 --------------------------------------------------------*/


#ifndef XSUBGRAPH_H
#define XSUBGRAPH_H


#include "vf2/argraph.h"

namespace vf2 {

Graph* ExtractSubgraph(Graph *g, node_id nodes, bool connected=true);

} // namespace vf2

#endif
