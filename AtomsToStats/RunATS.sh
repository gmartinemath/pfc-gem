#!/bin/bash
#This used to run RunATS.m on McGill's HPC Guillimin cluster
#PBS -l nodes=1:ppn=2
#PBS -l walltime=45:00:00
#PBS -V
#PBS -N GEM_ATS
#PBS -A its-555-aa
cd ~/GEM/AtomsToStats/
matlab -nosplash -nodesktop -nodisplay < ~/GEM/AtomsToStats/RunATS.m >> ~/GEM/AtomsToStats/GEM_PTS
