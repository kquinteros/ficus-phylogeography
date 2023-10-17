#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=144:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node
#SBATCH --mem=240G   # maximum memory per node
#SBATCH --job-name="phyluce_spades_ah"
#SBATCH --mail-user=kevinq@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluce167

phyluce_assembly_assemblo_spades \
  --conf /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/assembly_conf/assembly_Feb2020_ah.conf \
  --output /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/spades_assemblies \
  --log /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/logs \
  --cores 36
