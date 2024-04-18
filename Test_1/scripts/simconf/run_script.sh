#!/bin/sh

# SLURM options:

#SBATCH --job-name=Se82.2nubbsource_pads_bulk         	 # Job name
#SBATCH --partition=htc                  # Partition choice (most generally we work with htc, but for quick debugging you can use
										 #					 #SBATCH --partition=flash. This avoids waiting times, but is limited to 1hr)
#SBATCH --mem=16G                     	 # RAM
#SBATCH --licenses=sps                   # When working on sps, must declare license!!

#SBATCH --chdir=$TMPDIR
#SBATCH --time=2-0                 	 	 # Time for the job in format “minutes:seconds” or  “hours:minutes:seconds”, “days-hours”
#SBATCH --cpus-per-task=1                # Number of CPUs

# flsimulate
/sps/nemo/sw/snsw/2024/opt/falaise-5.1.1/bin/flsimulate -c /sps/nemo/scratch/mpetro/Projects/PhD/TrackReconstructionTests/Test_1/data/CAT/bb2nu/0/Se82.2nubb.conf -o SD.brio 

# flreconstruct reco pipeline 
/sps/nemo/sw/snsw/2024/opt/falaise-5.1.1/bin/flreconstruct -i SD.brio -p /pbs/home/m/mpetro/Projects/PhD/Codes/Job23_2nu_sensitivity/main/official-4.0.conf -o CD.brio

# flreconstruct SNCuts
/sps/nemo/sw/snsw/2024/opt/falaise-5.1.1/bin/flreconstruct -i CD.brio -p /pbs/home/m/mpetro/Projects/PhD/Codes/Job23_2nu_sensitivity/main/SNCutsPipeline.conf -o CDCut.brio

# flreconstruct MiModule
/sps/nemo/sw/snsw/2024/opt/falaise-5.1.1/bin/flreconstruct -p /pbs/home/m/mpetro/Projects/PhD/Codes/Job23_2nu_sensitivity/main/p_MiModule_v00.conf -i CDCut.brio 

# Default.root ---> data_file.root
# In this step I extract data from MiModule and save them in a new tree containing 
# the filtered events, their energies, angles, vertices, etc.
# This file is later analysed
/sps/nemo/sw/BxCppDev/opt/root-6.16.00/bin/root -l -b Job23.cpp