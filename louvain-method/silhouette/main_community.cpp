// File: main_community.cpp
// -- community detection, sample main file
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

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unistd.h>

#include "graph_binary.h"
#include "community.h"
#include "apsp.cpp"

using namespace std;

char   *filename_test = NULL;
char   *filename      = NULL;
char   *filename_w    = NULL;
char   *filename_part = NULL;
char   *filename_v    = NULL;
int    type           = UNWEIGHTED;
int    nb_pass        = 0;
double precision      = 0.000001;
int    display_level  = -2;
int    k1             = 16;
bool   verbose        = false;

void usage(char *prog_name, const char *more) {
	cerr << more;
	cerr << "usage: " << prog_name << " input_file [-w weight_file] [-p part_file] [-q epsilon] [-l display_level] [-v info_file] > tree_file" << endl << endl;
	cerr << "input_file: file containing the graph to decompose in communities." << endl;
	cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)." << endl;
	cerr << "-p file\tstart the computation with a given partition instead of the trivial partition." << endl;
	cerr << "\tfile must contain lines \"node community\"." << endl;
	cerr << "-q eps\ta given pass stops when the modularity is increased by less than epsilon." << endl;
	cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure." << endl;
	cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level." << endl;
	cerr << "-v file\tverbose mode: gives computation time, information about the hierarchy and modularity." << endl;
	exit(0);
}

void parse_args(int argc, char **argv) {
	if (argc<2)
		usage(argv[0], "Bad arguments number\n");

	for (int i = 1; i < argc; i++) {
		if(argv[i][0] == '-') {
			switch(argv[i][1]) {
				case 'w':
					type = WEIGHTED;
					filename_w = argv[i+1];
					i++;
					break;
				case 'p':
					filename_part = argv[i+1];
					i++;
					break;
				case 'q':
					precision = atof(argv[i+1]);
					i++;
					break;
				case 'l':
					display_level = atoi(argv[i+1]);
					i++;
					break;
				case 'k':
					k1 = atoi(argv[i+1]);
					i++;
					break;
				case 'v':
					verbose=true;
					filename_v = argv[i+1];
					i++;
					break;
				case 't':
					filename_test = argv[i+1];
					i++;
					break;
				default:
					usage(argv[0], "Unknown option\n");
			}
		} else {
			if (filename==NULL)
				filename = argv[i];
			else
				usage(argv[0], "More than one filename\n");
		}
	}
}

void display_time(const char *str, ofstream &foutput) {
	time_t rawtime;
	time ( &rawtime );
	foutput << str << ": " << ctime (&rawtime);
}

int** createMatrix(unsigned int nb_nodes, vector<float> weights, ofstream &foutput, char* filename) {
	int s, f, w;

	// initializing adjMatrix
	int** adjMatrix = new int*[nb_nodes];
	for (unsigned int i = 0; i < nb_nodes; ++i) {
		adjMatrix[i] = new int[nb_nodes];
	}
	for (int i = 0; i < nb_nodes; i++) { 
		for (int j = 0; j < nb_nodes; j++) {
			adjMatrix[i][j] = 0;
		}
	}

	ifstream fin;
	fin.open(filename,fstream::in);

	if (weights.size()!=0) { // weighted graph
		while (!fin.eof()) {
			fin >> s >> f >> w;
			adjMatrix[s][f] = w;
			adjMatrix[f][s] = w;
		}
	} else {
		while (!fin.eof()) {
			fin >> s >> f;
			adjMatrix[s][f] = 1;
			adjMatrix[f][s] = 1;
		}
	}

	fin.close();
	return adjMatrix;
}

int main(int argc, char **argv) {

	srand(time(NULL)+getpid());
	parse_args(argc, argv);

	ofstream foutput;
	foutput.open(filename_v, fstream::out);

	time_t time_begin, time_end;
	time(&time_begin);
	if (verbose)
		display_time("Begin", foutput);

	Community c(filename, filename_w, type, -1, precision);
	if (filename_part!=NULL) // start at a user-specified partition
		c.init_partition(filename_part); 

	foutput << "Created community " << endl; // TODO delete

	Graph g;
	bool improvement=true;
	double cov=c.coverage(), new_cov;
	int level=0;

	foutput << "Creating matrix" << endl; // TODO delete

	// create adjacency matrix
	int** adjMatrix = createMatrix(c.g.nb_nodes, c.g.weights, foutput, filename_test);

	// create apsp matrix
	int** apsp = floydWarshall(c.g.nb_nodes, adjMatrix);
	foutput << "----------------" << endl;
	for (int i = 0; i < c.g.nb_nodes; i++){ 
	  for (int j = 0; j < c.g.nb_nodes; j++) {
	    foutput << apsp[i][j] << " ";
	  }
	  foutput << endl;
	}
	foutput << "----------------" << endl;

	do {
		if (verbose) {
			foutput << "level " << level << ":\n";
			display_time("  start computation", foutput);
			foutput << "  network size: " 
							<< c.g.nb_nodes << " nodes, " 
							<< c.g.nb_links << " links, "
							<< c.g.total_weight << " weight." 
							<< endl;
		}

		// phase 1 - see community.h
		improvement = c.one_level();
		new_cov = c.coverage();

		if (++level==display_level) {
			g.display();
		}

		if (display_level==-1) {
			c.display_partition();
		}

		// phase 2
		g = c.partition2graph_binary();
		c = Community(g, -1, precision);

		if (verbose)
			foutput << "  coverage increased from " << cov << " to " << new_cov << endl;

		cov=new_cov;
		if (verbose)
			display_time("  end computation", foutput);

		if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
			improvement=true;
	} while(improvement);

	time(&time_end);
	if (verbose) {
		display_time("End", foutput);
		foutput << "Total duration: " << (time_end-time_begin) << " sec." << endl;
	}
	foutput << new_cov << endl;

	foutput.close();
}
