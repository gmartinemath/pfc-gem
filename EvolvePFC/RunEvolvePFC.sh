#!/bin/bash
#This used to run the RunEvolvePFC.m on McGill's HPC Guillimin cluster
#PBS -l nodes=1:ppn=12
#PBS -l walltime=50:00:00
#PBS -V
#PBS -N GEM_Evolve_PFC
#PBS -A its-555-aa
cd ~/GEM/EvolvePFC/
matlab -nosplash -nodesktop -nodisplay < ~/GEM/EvolvePFC/RunEvolvePFC.m >> ~/GEM/EvolvePFC/GEM_REPFC.out
