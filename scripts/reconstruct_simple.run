#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 100000
#SBATCH --cpus-per-task 5
#SBATCH --time 06:00:00
#SBATCH --exclusive

echo STARTING AT `date`

module purge
module load matlab
matlab -nosplash -nodisplay -nodesktop -r MAP2D_simple

cat /proc/cpuinfo | head -6

echo FINISHED at `date`

