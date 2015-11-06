// File: community-performance.cpp
// -- community detection using performance file
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
//            adapted 11/5/15 by Christina Tong
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time     : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "community-performance.h"

using namespace std;

// ----------------------------------------------------------------------------
// CONSTRUCTORS ---------------------------------------------------------------
Community::Community(char *filename, char *filename_w, int type, int nbp, double minp) {
  g = Graph(filename, filename_w, type); // create a binary graph
  size = g.nb_nodes;
  num_possible_edges = .5 * size * (size - 1);
  num_edges = nb_links / 2;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  n2c.resize(size); // vector; n2c[i] returns the community of node i
  in.resize(size); // vector; number of half-edges in a community i
  tot.resize(size); // vector; number of half-edges in or out of a community i

  // intialize all nodes into their own communities
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i; 
    tot[i] = g.weighted_degree(i);
    in[i]  = g.nb_selfloops(i);
  }

  nb_pass = nbp;
  min_perf = minp;
}

Community::Community(Graph gc, int nbp, double minp) {
  g = gc;
  size = g.nb_nodes;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  n2c.resize(size);
  in.resize(size);
  tot.resize(size);

  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    tot[i] = g.weighted_degree(i);
    in[i]  = g.nb_selfloops(i);
  }

  nb_pass = nbp;
  min_modularity = minp;
}

// ----------------------------------------------------------------------------
// FUNCTION DEFINITIONS -------------------------------------------------------

void Community::init_partition(char * filename) {
  ifstream finput;
  finput.open(filename,fstream::in);

  // read partition
  while (!finput.eof()) {
    unsigned int node, comm;
    finput >> node >> comm;
    
    if (finput) {
      int old_comm = n2c[node];
      neigh_comm(node);

      remove(node, old_comm, neigh_weight[old_comm]);

      unsigned int i=0;
      for ( i=0 ; i<neigh_last ; i++) {
        unsigned int best_comm = neigh_pos[i];
        float best_nblinks = neigh_weight[neigh_pos[i]];

        if (best_comm==comm) {
          insert(node, best_comm, best_nblinks);
          break;
        }
      }

      if (i==neigh_last)
        insert(node, comm, 0);
    }
  }
  finput.close();
}

// print each node, its community, 
// the number of half-edges in the community,
// and the number of half-edges in or out of the community
void Community::display() {
  for (int i=0 ; i<size ; i++)
    cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i] ;
  
  cerr << endl;
}

/* 
* @return - the performance score of the current graph partition
*/
double Community::performance() {
    double perf = 0;

  /* MY PSEUDOCODE:

  for every vertex u in the graph 
    for every vertex v > u (strict inequality to prevent duplicates)
        if u and v are in the same community
          if the edge (u,v) exists,
            f = f + w(u,v)
          else do nothing
        else 
          if the edge between them does not exist, 
            g = g + w (where w is 1 for unweighted graph, 
                      otherwise it's the average weight for nonexistent edges)
          else do nothing

  add f and g
  divide by num_possible_edges
  */

    return perf;
}

/* 
* Compute the performance of the initial partition of G,
* in which every node is its own cluster. This is a special 
* case of calculating performance. That is, f(G) = 0 
* because every node is one clustering and thus there are 
* no internal edges within communities. 
*
* g(G) = num_possible_edges - |E|, because at this point
* every nonexistent edge (disregarding loops) between 
* 2 vertices is between two separate communities, so they
* count in the g(G) score.
* 
* If this a weighted graph, we multiply the number of nonexistent 
* external edges by the average weight of the nonexistent external edges,
* which we calculate from the other weighted edges.
*
* Finally, we divide g(G) by the total number of possible edges
* between all vertices (num_possible_edges).
*
* @return perf - performance score
*/
double Community::initial_performance() {

    int g = num_possible_edges - num_edges;   // TODO: how does this work for later iterations of the algorithm?

/* Alternate method to find the number of 
* nonexistent external edges (g(G) for unweighted graph): 
* count the number of 0's strictly above the diagonal 
* in the adjacency matrix A.

for (int i = 0; i < size; i++) {
    for (int j = i+1; i < size; j++) {
        if (A[i][j] == 0) g++;
    }
}
*/

    // TODO: calculate avg weight, for weighted graphs, here

    return ((float)(g)) / ((float)(num_possible_edges));
}

// compute the neighboring communities of the node
void Community::neigh_comm(unsigned int node) {
  for (unsigned int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1; // set the weight of each node in 
  
  neigh_last=0;
  pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(node);
  unsigned int deg = g.nb_neighbors(node);

  neigh_pos[0]=n2c[node];
  neigh_weight[neigh_pos[0]]=0;
  neigh_last=1;

  for (unsigned int i=0 ; i<deg ; i++) {
    unsigned int neigh        = *(p.first+i);
    unsigned int neigh_comm   = n2c[neigh];
    double neigh_w = (g.weights.size()==0)?1.:*(p.second+i);
    
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
        neigh_weight[neigh_comm]=0.;
        neigh_pos[neigh_last++]=neigh_comm;
      }
      neigh_weight[neigh_comm]+=neigh_w;
    }
  }
}

void Community::partition2graph() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  for (int i=0 ; i<size ; i++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(i);

    int deg = g.nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
    }
  }
}

// print the graph partition
void Community::display_partition() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  for (int i=0 ; i<size ; i++)
    cout << i << " " << renumber[n2c[i]] << endl;
}


Graph Community::partition2graph_binary() {
  // Renumber communities
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  // Compute communities
  vector<vector<int> > comm_nodes(final);
  for (int node=0 ; node<size ; node++) {
    comm_nodes[renumber[n2c[node]]].push_back(node);
  }

  // Compute weighted graph
  Graph g2;
  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(comm_nodes.size());

  int comm_deg = comm_nodes.size();
  for (int comm=0 ; comm<comm_deg ; comm++) {
    map<int,float> m;
    map<int,float>::iterator it;

    int comm_size = comm_nodes[comm].size();
    for (int node=0 ; node<comm_size ; node++) {
      pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(comm_nodes[comm][node]);
      int deg = g.nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
        int neigh        = *(p.first+i);
	      int neigh_comm   = renumber[n2c[neigh]];
	      double neigh_weight = (g.weights.size()==0)?1.:*(p.second+i);

	      it = m.find(neigh_comm);
	      if (it==m.end())
          m.insert(make_pair(neigh_comm, neigh_weight));
        else
          it->second+=neigh_weight;
      }
    }
    g2.degrees[comm]=(comm==0)?m.size():g2.degrees[comm-1]+m.size();
    g2.nb_links+=m.size();

    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight  += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
  }

  return g2;
}

// Computes community partition for one pass
bool Community::one_level() {
  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  //double new_mod   = modularity();
  //double cur_mod   = new_mod;

  // chose the random order in which we examine the vertices
  vector<int> random_order(size);
  for (int i=0 ; i<size ; i++)
    random_order[i]=i;
  for (int i=0 ; i<size-1 ; i++) {
    int rand_pos = rand()%(size-i)+i;
    int tmp      = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // potential increase from placing a node into a neighboring community
  double increase;

  /*
  * repeat while 
  *   there is an improvement of performance
  *   or there is an improvement of performance greater than a given epsilon 
  *   or a predefined number of passes have been completed
  */
  do {
    // cur_mod = new_mod;
    nb_moves = 0;   // number of nodes moved into new communities
    nb_pass_done++;   // count the number of passes

    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
        // int node = node_tmp;
        int node = random_order[node_tmp];
        int node_comm     = n2c[node];
        double w_degree = g.weighted_degree(node);

        // computation of all neighboring communities of current node
        neigh_comm(node);
        
        // remove node from its current community
        remove(node, node_comm, neigh_weight[node_comm]);
        
        // calculate the change in the numerator due to removal of node
        remove_perf_change = remove_perf_change(int node, int comm, double dnodecomm);

        // compute the nearest community for node
        // default choice for future insertion is the former community 
        int best_comm        = node_comm;
        double best_nblinks  = 0.;
        double best_increase = 0.;

        // consider putting the node into each of its neighboring communities
        for (unsigned int i=0 ; i<neigh_last ; i++) {
          increase = performance_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);


          // TODO: deal with the case where it's best to just leave the node in its own community?
          // if a new best increase is calculated, update accordingly
          if (increase>best_increase) {
            best_comm     = neigh_pos[i];
            best_nblinks  = neigh_weight[neigh_pos[i]];
            best_increase = increase;
          }
        }

        // insert node in the nearest community
        insert(node, best_comm, best_nblinks);
       
        // save old performance for calculation in while loop boolean guard
        old_perf = cur_perf;

        // update performance accordingly; change will be >= 0
        // remove_perf_change can be negative
        cur_perf = cur_perf + remove_perf_change + best_increase;

        // if we put the node in a different community, we made an improvement
        if (best_comm!=node_comm) {
          nb_moves++;
        }
    }

    // TODO what is this doing?
    double total_tot=0;
    double total_in=0;
    for (unsigned int i=0 ; i<tot.size() ;i++) {
      total_tot+=tot[i];
      total_in+=in[i];
    }

    if (nb_moves>0)
      improvement=true;
    
  } while (nb_moves>0 && cur_perf-old_perf>min_perf); // TODO: change

  return improvement;
}

