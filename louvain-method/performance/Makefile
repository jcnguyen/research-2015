#!/bin/bash
CC=g++
CFLAGS=-ansi -O3 -Wall
LDFLAGS=-ansi -lm -Wall
EXEC=convert community hierarchy display
OBJ1=graph.o
OBJ2=graph_binary.o community.o

all: $(EXEC)

convert : $(OBJ1) main_convert.o
	$(CC) -o $@ $^ $(LDFLAGS)

community : $(OBJ2) main_community.o
	$(CC) -o $@ $^ $(LDFLAGS)

hierarchy : main_hierarchy.o
	$(CC) -o $@ $^ $(LDFLAGS)

display : main_display.o
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
