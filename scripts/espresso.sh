#!/bin/bash

#SBATCH --job-name=anytime-pmnk-landscapes

##SBATCH --mail-user=
##SBATCH --mail-type=ALL

##SBATCH --nodes=1
##SBATCH --partition normal
##SBATCH --qos=general-compute

##SBATCH --mem-per-cpu=
##SBATCH --mem=

##SBATCH --time=0-00:05:00
##SBATCH --ntasks-per-node=12
##SBATCH --ntasks=1
##SBATCH --requeue

srun -o "$1" -e "$2" "$3" "${@:4:$#-1}" --oversubscribe -t=05:00

exit 0
