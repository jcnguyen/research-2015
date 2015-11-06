// File: community-performance.h
// -- community detection header file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large 
// networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above 
// mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef COMMUNITY_PERFORMANCE_H
#define COMMUNITY_PERFORMANCE_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

#include "graph_binary.h"

using namespace std;

class Community {
 public:  

  // --------------------------------------------------------------------------
  // GLOBAL VARIABLES ---------------------------------------------------------

  // TODO unsure if jennifer's comments correct?
  vector<double> neigh_weight; // weight of the neighboring communities all of vertices in this community; index = vertex, double = weight 
  vector<unsigned int> neigh_pos; // the current community that a vertex is in; index = vertex, neigh_pos[i] = communityID
  unsigned int neigh_last; // the position of the final neighbor of this community

  Graph g;            // network to compute communities for
  int size;           // number of nodes in the network and size of all vectors
  int num_possible_edges; // total possible number of edges between all vertices
  int num_edges;    // the number of actual edges in g
  int remove_perf_change;     // change in performance due to removal of node from a community
  vector<int> n2c;    // community to which each node belongs

  // used to compute the modularity participation of each community
  vector<double> in;  // number of half-edges in a community i
  vector<double> tot; // number of half-edges in or out a community i

  // number of pass for one level computation
  // if -1, compute as many pass as needed to increase modularity
  int nb_pass;

  // a new pass is computed if the last one generated 
  // an increase greater than min_perf
  // if min_perf = 0, a minor increase is enough to go for one more pass
  double min_perf;

  // --------------------------------------------------------------------------
  // CONSTRUCTORS -------------------------------------------------------------

  // reads graph from file using graph constructor
  // type defines the weighted/unweighted status of the graph file
  // nb_pass is the maximum number of passes of the algorithms
  Community (char *filename, char *filename_w, int type, int nb_pass, double min_perf);
  Community (Graph g, int nb_pass, double min_perf); // copy graph

  // --------------------------------------------------------------------------
  // FUNCTION DECLARATIONS ----------------------------------------------------

  // initiliazes the partition based on a file specification rather than with
  // the default of each node as its own community
  void init_partition(char *filename_part);

  // display the community of each node
  void display();

  // remove node from its current community with which it has dnodecomm links
  inline void remove(int node, int comm, double dnodecomm);

  // insert node in the community with which it shares dnodecomm links
  inline void insert(int node, int comm, double dnodecomm);

  inline void Community::remove_perf_change(int node, int comm, double dnodecomm);

  // compute the gain of performance if node was inserted in the community,
  // given that the node has dnodecomm links to the community. The formula is:
  // TODO: formula
  // where 
  //  In(comm)    = number of half-links strictly inside comm
  //  Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
  //  d(node,com) = number of links from node to comm
  //  deg(node)   = node degree
  //  m           = number of links (edges)
  inline double performance_gain(int node, int comm, double dnodecomm, double w_degree);

  // compute the set of neighboring communities of a node
  // for each community, gives the number of links from node to comm
  void neigh_comm(unsigned int node);

  // compute the performance of the current graph partition
  double performance();

  // displays the graph of communities as computed by one_level
  void partition2graph();

  // displays the current partition (with communities renumbered from 0 to k-1)
  void display_partition();

  // generates the binary graph of communities as computed by one_level
  Graph partition2graph_binary();

  // compute communities of the graph for one level
  // return true if some nodes have been moved
  bool one_level();
};

// ----------------------------------------------------------------------------
// INLINE FUNCTION DEFINITIONS ------------------------------------------------

inline void Community::remove(int node, int comm, double dnodecomm) {
  assert(node>=0 && node<size);

  tot[comm] -= g.weighted_degree(node); // subtract from # of external edges
  in[comm]  -= 2*dnodecomm + g.nb_selfloops(node); // subtract from # of internal edges
  n2c[node]  = -1; // not in a community
}

/*
* @param node - the node to be removed
* @param comm - community containing node
* @param dnodecomm - number of links from node to comm
* 
* @return performance change from removing node from comm
*/
inline void Community::remove_perf_change(int node, int comm, double dnodecomm) {
  assert(node>=0 && node<size);
  int change = 0;
  
  // look at all possible edges, existent or nonexistent, with one end at node
  for(int i=0; i<size; i++) {

      // Removing node from comm only changes its relations with
      // vertices in comm. To all other vertices, it is still just
      // some vertex in a different community.
      if (n2c[i] == comm) {

          // since node and vertices v in comm are no longer in the same
          // community, we increase g(G) for nonexistent edges (node,v),
          if (A[node][i] == NO_EDGE) {
              change += 1;    // TODO: change for weighted graphs
          } 

          // removing node from its community deducts from f(G)
          // those edges that went from (node, v), where v is a vertex in comm 
          else {
              change -= A[node][i];
          }

      }
  }

  return (float)change/(float)num_possible_edges;
}

/* 
* computes the potential performance gain
*
* TODO: remove some of these inputs?
* @param node - node to place
* @param comm - community to put node into
* @param dnodecomm - number of links from node to comm
* @param w_degree - node degree
* @return - change in performance 
*/
inline double Community::performance_gain(int node, int comm, double dnodecomm, double w_degree) {
  assert(node>=0 && node<size);
  int change = 0;
  
  // look at all possible edges, existent or nonexistent, with one end at node
  for(int i=0; i<size; i++) {

      // Adding a node into comm only changes its relations with
      // vertices in comm. To all other vertices, it is still just
      // some vertex in a different community.
      if (n2c[i] == comm) {

          // subtract from g(G) the nonexistent edges (node, v),
          // where v is in comm. We don't count nonexistent edges
          // within the community toward our performance score
          if (A[node][i] == NO_EDGE) {
              change -= 1;    // TODO: change for weighted graphs
          } 

          // adding node to comm adds to f(G)
          // for existent edges (node, v) where v is a vertex in comm
          else {
              change += A[node][i];
          }
      }
  }

  return (float)change/(float)num_possible_edges;
}

#endif 
// COMMUNITY_PERFORMANCE_H
