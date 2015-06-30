/* 
HOW TO USE
----------
compiling:
make

running:
unweighted: ./indexzero -i input_filename.txt > output_filename.txt
weighted:   ./indexzero -i input_filename.txt -w > output_filename.txt

*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

char *infile     = NULL; 
int  type        = UNWEIGHTED;

void convert_to_index_0(char *filename, int type) {
  // opens up file reader
  ifstream finput;
  finput.open(filename,fstream::in);

  unsigned int src, dest; // starting vertex and ending vertex
  unsigned int new_src, new_dest;
  double weight=1.; // weight between src and dest vertex

  if (type == WEIGHTED) {
    while(finput >> src >> dest >> weight) {
      new_src = src - 1;
      new_dest = dest - 1;
      cout << new_src << " " << new_dest << " " << weight << endl;
    }
  } else if (type == UNWEIGHTED){
    while(finput >> src >> dest) {
      new_src = src - 1;
      new_dest = dest - 1;
      cout << new_src << " " << new_dest << endl;
    }
  }

  finput.close();
}

void usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " -i input_file -o outfile [-r] [-w outfile_weight] [-d]" << endl << endl;
  cerr << "read the graph and convert it to binary format." << endl;
  cerr << "-r\tnodes are renumbered from 0 to nb_nodes-1 (the order is kept)." << endl;
  cerr << "-w\tread the graph as a weighted one and writes the weights in a separate file." << endl;
  cerr << "-d\tdisplays graph." << endl;
  exit(0);
}

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
        case 'w' :
          type = WEIGHTED;
          i++;
          break;
        default:
          usage(argv[0], "Unknown option\n");
      }
    } else {
      usage(argv[0], "More than one filename\n");
    }
  }
  if (infile==NULL)
    usage(argv[0], "Infile missing\n");
}

// Uses txt file to build a graph
int main(int argc, char **argv) {
  parse_args(argc, argv);

  convert_to_index_0(infile, type);
}

