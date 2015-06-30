// File: main_convert.cpp
// -- conversion of a graph from ascii to binary, sample main file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time     : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "graph.h"

using namespace std;

char *infile     = NULL; // file of ascii graph
char *outfile    = NULL; // file of binary graph
char *outfile_w  = NULL; // file of weights
int  type        = UNWEIGHTED; // graph is unweighted by default
bool do_renumber = false; // renumber the vertices? 
bool display     = false; // display graph?

// How to use this exe
void usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " -i input_file -o outfile [-r] [-w outfile_weight] [-d]" << endl << endl;
  cerr << "read the graph and convert it to binary format." << endl;
  cerr << "-r\tnodes are renumbered from 0 to nb_nodes-1 (the order is kept)." << endl;
  cerr << "-w filename\tread the graph as a weighted one and writes the weights in a separate file." << endl;
  cerr << "-d\tdisplays graph." << endl;
  exit(0);
}

// Parses command-line arguments
void parse_args(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
        case 'i':
          if (i==argc-1)
            usage(argv[0], "Infile missing\n");
          infile = argv[i+1];
          i++;
          break;
        case 'o':
          if (i==argc-1)
            usage(argv[0], "Outfile missing\n");
          outfile = argv[i+1];
          i++;
          break;
        case 'w' :
          type = WEIGHTED;
          outfile_w = argv[i+1];
          i++;
          break;
        case 'r' :
          do_renumber = true;
          break;
        case 'd' :
          display = true;
          break;
        default:
          usage(argv[0], "Unknown option\n");
      }
    } else {
      usage(argv[0], "More than one filename\n");
    }
  }
  if (infile==NULL || outfile==NULL)
    usage(argv[0], "In or outfile missing\n");
}

// Uses txt file to build a graph
int main(int argc, char **argv) {
  parse_args(argc, argv);

  Graph g(infile, type);
  g.clean(type);

  if (do_renumber)
    g.renumber(type);

  if (display)
    g.display(type);
  
  g.display_binary(outfile, outfile_w, type);

}
