#!/bin/sh

numactl -H | grep "nodes" | awk '{print $2}'
