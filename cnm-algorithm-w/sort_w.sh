#!/bin/bash

NAME_LIST=(
	"lesmis-cnm-modularity.groups"
	"football-cnm-modularity.groups"
	"celegansneural-cnm-modularity.groups"
)

for name in "${NAME_LIST[@]}"
do
	sort -n < "${name}.groups" > "${name}.graph"
done