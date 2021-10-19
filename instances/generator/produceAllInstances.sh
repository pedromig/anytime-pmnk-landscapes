#!/bin/bash

# the directory to put the instances
dirInstances=".."

# the directory where the generator can be found
dirGenerator=.

# list of instance parameters to produce
# of course you can change this list
listRho="-0.7 -0.3 0.0 0.3 0.7"
listM="2 3 5 7"
listN="18 32 64 128"
listK="1 2 4 8"

nbInst=5

# create the directory of instances if necessary

if [ ! -d ${dirInstances} ]; then mkdir ${dirInstances}; fi

# Go to produce all the instances.
# It can take a lot of time (hours) if there is a lot of large instances !
for M in $listM; do
	for N in $listN; do
		for K in $listK; do
			if [ "$(echo "$K" '<' "$N" | bc -l)" -eq 1 ]; then
				for r in $listRho; do
					seed=0
					for ((n = 0; n < ${nbInst}; n++)); do
						name=rmnk_${r}_${M}_${N}_${K}_${seed}.dat
						if [ "$(echo "$r" '> (-1.0 / ('"$M"' - 1))' | bc -l)" -eq 1 ]; then
							R --slave --no-restore --file=${dirGenerator}/rmnkGenerator.R --args $r $M $N $K $seed ${dirInstances}/${name}
							# you can use this command if you need: ${dirGenerator}/rmnkGenerator.R $r $M $N $K $n ${dirInstances}/${name}
						fi
						seed=$((seed + 1))
					done
				done
			fi
		done
	done
done
