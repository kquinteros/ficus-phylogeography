#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=24   # 24 processor core(s) per node
#SBATCH --mem=50G   # maximum memory per node
#SBATCH --job-name="UCE_petiolaris "
#SBATCH --mail-user=kevinq@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluce167


illumiprocessor \
    --input raw-fastq/ \
    --output clean-fastq \
    --config illumiprocessor.conf \
    --cores 4

#######################################
###########Quality control#############
#######################################

# move to the directory holding our cleaned reads
cd clean-fastq/

# run this script against all directories of reads

for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
done

#make sure we are at the top-level of our uce tutorial directory
cd uce-tutorial


module unuse /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module use /opt/rit/modules

module load java/1.7.0_55
module load bowtie/1.1.2

# run the assembly
phyluce_assembly_assemblo_trinity \
    --conf assembly.conf \
    --output trinity-assemblies \
    --clean \
    --cores 12
