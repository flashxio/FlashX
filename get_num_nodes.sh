#!/bin/sh

num_nodes=`numactl -H | grep "nodes" | awk '{print $2}'`
if [ "$num_nodes" == 0 ]; then
	echo 1
elif [ -z $num_nodes ]; then
	echo 1
else
	echo $num_nodes
fi
