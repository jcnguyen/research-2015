 // File: community-coverage.cpp
// -- community detection using coverage file
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

#include "community.h"

using namespace std;

// ----------------------------------------------------------------------------
// CONSTRUCTORS ---------------------------------------------------------------
Community::Community(char *filename, char *filename_w, int type, int nbp, double minc) {
  g = Graph(filename, filename_w, type); // create a graph
  size = g.nb_nodes;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  n2c.resize(size);
  in.resize(size);
  tot.resize(size);

  for (int i=0 ; i<size ; i++) {
    n2c[i] = i; // initially, each node is its own community
    tot[i] = g.weighted_degree(i);
    in[i]  = g.nb_selfloops(i);
  }

  nb_pass = nbp;
  min_coverage = minc;
}

Community::Community(Graph gc, int nbp, double minc) {
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
  min_coverage = minc;
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

void Community::display() {
  for (int i=0 ; i<size ; i++)
    cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i] ;
  
  cerr << endl;
}

double Community::silhouetteAux(int** apsp) {
  // initialize community matrix
  int** whatCommunity = new int*[size];
  for (int i = 0; i < size; ++i) {
    whatCommunity[i] = new int[n2c.size];
  }

  // add vertex into corresponding community vector
  for (i = 0; i < n2c.size(); i++) {
    whatCommunity[n2c[i]].push_back(i);
  }
  //--------------------------------------------------------------------//

  vector<double> allCommunitySilhouettes[size];

  // for every community
  for (c = 0; c < size; c++) {
    sizeOfCommunity = whatCommunity[c].size();
    vector<double> silhouetteWidths;

    // for every node in community c
    for (i = 0; i < sizeOfCommunity; i++) {
      double averageInnerDistance = 0.0;

      // find avg of shortest distance between v_i and every other vertex in same community
      if (sizeOfCommunity == 1) { // singleton
        averageInnerDistance = 1.0;
      } else {
        for (k = 0; k < sizeOfCommunity; k++) {
          vertexK = whatCommunity[c][k];
          if (i != vertexK) {
            averageInnerDistance += apsp[i][vertexK]
          }
        }
        averageInnerDistance = averageInnerDistance/(sizeOfCommunity - 1);
      }

      vector<double> distanceOptions[sizeOfCommunity];
      distanceOptions[c] = 0;

      // finds distance between v_i and every vertex in different community
      for (n = 0; n < size; n++) {
        double currentDistance = 0;

        // not in our community c
        if (n != c) {
          sizeOfCommunityN =  whatCommunity[n].size();

          for (k = 0; k < sizeOfCommunityN; k++) {
            vertexK = whatCommunity[n][k];
            currentDistance += apsp[i][vertexK];
          }
          currentDistance = currentDistance/(sizeOfCommunityN - 1)
          distanceOptions[n] = currentDistance;
        }

      }

      double minimumOfAverageOutsideDistance = min_element( distanceOptions.begin(), distanceOptions.end() );
      double maxOfDistances = 0;

      if {averageInnerDistance > minimumOfAverageOutsideDistance}
        maxOfDistances = averageInnerDistance;
      else
        maxOfDistances = minimumOfAverageOutsideDistance;

      silhouetteWidth[i] = (minimumOfAverageOutsideDistance - averageInnerDistance)/maxOfDistances; 

    }

    // calculate the silhouette width for the entrie community
    double communitySilhouette = 0;
    for (k = 0; k < sizeOfCommunity; k++) {
      communitySilhouette += silhouetteWidth[k]
    }

    allCommunitySilhouettes[c] = communitySilhouette/sizeOfCommunity

  }

  double globalSilhouetteIndex = 0;
  for (k = 0; k < size; k++) {
    globalSilhouetteIndex += allCommunitySilhouettes[k];
  }

  return globalSilhouetteIndex/size
}

double Community::coverage() {
  double cov  = 0.;
  double tot_weight = (double)g.total_weight;

  for (int i=0 ; i<size ; i++) {
    if (tot[i]>0)
      cov += (double)in[i]/tot_weight;
  }

  return cov;
}

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
  bool   improvement  = false;
  int    nb_moves;
  int    nb_pass_done = 0;
  double new_cov      = coverage();
  double cur_cov      = new_cov;

  vector<int> random_order(size);
  for (int i=0 ; i<size ; i++)
    random_order[i]=i;

  for (int i=0 ; i<size-1 ; i++) {
    int rand_pos = rand()%(size-i)+i;
    int tmp      = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // repeat while 
  //   there is an improvement of coverage
  //   or there is an improvement of coverage greater than a given epsilon 
  //   or a predefined number of pass have been done
  do {
    cur_cov = new_cov;
    nb_moves = 0;
    nb_pass_done++;

    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
//      int node = node_tmp;
      int node = random_order[node_tmp];
      int node_comm     = n2c[node];
      double w_degree = g.weighted_degree(node);

      // computation of all neighboring communities of current node
      neigh_comm(node);
      // remove node from its current community
      remove(node, node_comm, neigh_weight[node_comm]);

      // compute the nearest community for node
      // default choice for future insertion is the former community (TODO is this tiebreaking?)
      int best_comm        = node_comm;
      double best_nblinks  = 0.;
      double best_increase = 0.;
      for (unsigned int i=0 ; i<neigh_last ; i++) {
        double increase = coverage_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
        // cerr << "increase: " << increase << ", best: " << best_increase<< endl; // todo
        if (increase>=best_increase) {
          best_comm     = neigh_pos[i];
          best_nblinks  = neigh_weight[neigh_pos[i]];
          best_increase = increase;
        }
      }

      // insert node in the nearest community
      insert(node, best_comm, best_nblinks);
     
      if (best_comm!=node_comm)
        nb_moves++;
    }

    double total_tot=0;
    double total_in=0;
    for (unsigned int i=0 ; i<tot.size() ;i++) {
      total_tot+=tot[i];
      total_in+=in[i];
    }

    new_cov = coverage();
    if (nb_moves>0)
      improvement=true;
    
  } while (nb_moves>0 && new_cov-cur_cov>min_coverage);

  return improvement;
}
