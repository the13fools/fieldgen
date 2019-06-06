#!/usr/bin/env bash

for i in ../fieldbenchmark/RFGParametrizationBenchmark/*; do
	for j in $(seq 0.00 0.2 1.00)
	do
    ./fieldgen $i $i --degree=4 --alignToCurvature --bisectT --sampleToFaces --s=1 --t=$j
	done

  echo $i
done

# ./fieldgen data/rocker_arm.obj data/rocker_arm.repv --degree=4 --alignToCurvature --bisectT --s=1 --t=.75
