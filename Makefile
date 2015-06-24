#!/bin/bash

CC      = g++
CFLAGS  = -ansi -O5 -Wall
LDFLAGS = -ansi -lm -Wall
EXEC    = convert community_modularity community_coverage hierarchy indexzero
OBJ1    = graph.o
OBJ2    = graph_binary.o community_modularity.o
OBJ3    = graph_binary.o community_coverage.o


all: $(EXEC)

convert : $(OBJ1) main_convert.o
	$(CC) -o $@ $^ $(LDFLAGS)

community_modularity : $(OBJ2) main_community_modularity.o
	$(CC) -o $@ $^ $(LDFLAGS)

community_coverage : $(OBJ3) main_community_coverage.o
	$(CC) -o $@ $^ $(LDFLAGS)

hierarchy : main_hierarchy.o
	$(CC) -o $@ $^ $(LDFLAGS)

indexzero : main_indexzero.o
	$(CC) -o $@ $^ $(LDFLAGS)

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o *~ $(EXEC)
