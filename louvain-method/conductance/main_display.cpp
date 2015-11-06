#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

char *filename = NULL;

// How to use this exe
void usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file > graph_file" << endl << endl;
  cerr << "input_file: read the community tree from this file." << endl;
  cerr << "graph_file: store the final community graph in this file." << endl;
  exit(0);
}

// Parses arguments
void parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if (filename==NULL)
      filename = argv[i];
    else
      usage(argv[0], "More than one filename\n");
  }
  
  if (filename==NULL)
    usage(argv[0], "No input file has been provided.\n");
}

// Main sequence
int main(int argc, char **argv) {
  parse_args(argc, argv);

  // lists the communities in which the node belongs to at that particular levell (iteration),
  // where index = level, vector<int>[level] = node id, vector<vector><int>[level][node] = community
  // note that vector<int> will be smaller as iterations pass
  vector<vector<int> >levels; 

  ifstream finput;
  finput.open(filename,fstream::in);

  // construct the levels
  int l=-1;
  while (!finput.eof()) {
    int node, community;
    finput >> node >> community;

    if (finput) {
      if (node==0) {
        l++;
        levels.resize(l+1);
      }
      levels[l].push_back(community);
    }
  }

  // prints the info
  // note that levels[0] is when each node is its own community,
  // so levels[0].size is the total number of nodes
  vector<int> n2c(levels[0].size());

  // each node (index) is its own community
  for (unsigned int node=0 ; node<levels[0].size() ; node++)
    n2c[node]=node;
    
  // setting the community in which the node belongs
  for (l=0 ; l<(levels.size() - 1) ; l++)
    for (unsigned int node=0 ; node<levels[0].size() ; node++)
      n2c[node] = levels[l][n2c[node]];
   
  // print "[node] [community]"
  for (unsigned int node=0 ; node<levels[0].size() ; node++)
    cout << node << " " << n2c[node] << endl;
}
