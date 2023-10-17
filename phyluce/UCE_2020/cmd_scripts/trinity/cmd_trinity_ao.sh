#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=216:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mem=50G   # maximum memory per node
#SBATCH --job-name="phyluce_trinity_ao"
#SBATCH --mail-user=kevinq@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluce162
​
module unuse /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module use /opt/rit/modules
​
module load java/1.7.0_55
module load bowtie/1.1.2
​
​
phyluce_assembly_assemblo_trinity \
    --conf /ptmp/kevinq/assembly_conf/assembly_Feb2020_ao.conf \
    --output /ptmp/kevinq/trinity-assemblies \
    --log /ptmp/kevinq/logs \
    --clean \
    --cores 16
