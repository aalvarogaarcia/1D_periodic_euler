#!/bin/bash

for n in $(seq 3 30); do
	N=$((2**$n))
	echo $N $(taskset -c 0 ./comparisonVector.p $N)
done