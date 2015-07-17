// File: community-coverage.h
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
// Time     : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef COMMUNITY_COVERAGE_H
#define COMMUNITY_COVERAGE_H

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

  // TODO wtf are these; fix comments
  vector<double> neigh_weight; // weight of the neighboring communities all of vertices in this community; index = vertex, double = weight 
  vector<unsigned int> neigh_pos; // the current community that a vertex is in; index = vertex, neigh_pos[i] = communityID
  unsigned int neigh_last; // the position of the final neighbor of this community

  Graph g;            // network to compute communities for
  int size;           // number of nodes in the network and size of all vectors
  vector<int> n2c;    // community to which each node belongs

  // used to compute the modularity participation of each community
  vector<double> in;  // number of half-edges in a community i
  vector<double> tot; // number of half-edges in or out a community i

  // number of pass for one level computation
  // if -1, compute as many pass as needed to increase coverage
  int nb_pass;

  // a new pass is computed if the last one has generated 
  // an increase greater than min_coverage
  // if min_coverage = 0, a minor increase is enough to go for one more pass
  double min_coverage;

  // --------------------------------------------------------------------------
  // CONSTRUCTORS -------------------------------------------------------------

  // reads graph from file using graph constructor
  // type defined the weighted/unweighted status of the graph file
  Community (char *filename, char *filename_w, int type, int nb_pass, double min_coverage);
  Community (Graph g, int nb_pass, double min_coverage); // copy graph

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

  // compute the gain of coverage if node was inserted in the community,
  // given that the node has dnodecomm links to the community. The formula is:
  inline double coverage_gain(int node, int comm, double dnodecomm, double w_degree);

  // compute the set of neighboring communities of node
  // for each community, gives the number of links from node to comm
  void neigh_comm(unsigned int node);

  // compute the coverage of the current graph partition
  double coverage();

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

  tot[comm] -= g.weighted_degree(node);
  in[comm]  -= 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]  = -1;
}

inline void Community::insert(int node, int comm, double dnodecomm) {
  assert(node>=0 && node<size);

  tot[comm] += g.weighted_degree(node);
  in[comm]  += 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]=comm;
}

inline double Community::coverage_gain(int node, int comm, double dnodecomm, double w_degree) {
  // assert(node>=0 && node<size);
  // double tot_weight = (double)g.total_weight;
  // return inc/tot_weight;

  assert(node>=0 && node<size);

  double inc = (double)in[comm];
  double degc = (double)w_degree;
  double m2   = (double)g.total_weight;
  double dnc  = (double)dnodecomm;
  
  return (dnc - inc*degc/m2);
}


#endif // COMMUNITY_COVERAGE_H
