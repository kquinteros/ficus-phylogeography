#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH -t 120:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "phyluce_trinity_ab"   # job name
#SBATCH --mail-user=jsatler@iastate.edu   # email address
#SBATCH --mail-type=END



# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluce

module unuse /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module use /opt/rit/modules

module load java/1.7.0_55
module load bowtie/1.1.2


phyluce_assembly_assemblo_trinity \
    --conf /ptmp/LAS/phylo-lab/jsatler/phyluce/assembly_conf/assembly_rd23_ab.conf \
    --output /ptmp/LAS/phylo-lab/jsatler/phyluce/trinity-assemblies \
    --clean \
    --cores 16
