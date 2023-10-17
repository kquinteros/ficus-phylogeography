#!/bin/bash
# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=144:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node
#SBATCH --mem=240G   # maximum memory per node
#SBATCH --job-name="clean_align"
#SBATCH --mail-user=kevinq@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluceBio

cd /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE

# run this script against all directories of reads

for i in spades_assemblies/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

################################################################################
                           #Finding UCE Loci
################################################################################

phyluce_assembly_match_contigs_to_probes \
    --contigs /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/spades_assemblies/contigs \
    --probes /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/hymenoptera-v2-PRINCIPAL-bait-set.fasta.txt \
    --output /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/uce-search-results \
    --verbosity INFO
################################################################################
                          #Extracting UCE Loci
################################################################################
# create an output directory for this taxon set - this just keeps
# things neat

#mkdir -p taxon-sets/all

# create the data matrix configuration file
phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf

# change to the taxon-sets/all directory
cd taxon-sets/all

# make a log directory to hold our log files - this keeps things neat
#mkdir log

# get FASTA data for taxa in our taxon set
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../spades_assemblies/contigs \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log

################################################################################
                      #Exploding the monolithic FASTA File
################################################################################

# explode the monolithic FASTA by taxon (you can also do by locus)
phyluce_assembly_explode_get_fastas_file \
    --alignments all-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon

# get summary stats on the FASTAS
for i in exploded-fastas/*.fasta;
do
  phyluce_assembly_get_fasta_lengths --input $i --csv >> exploded.txt;
done

################################################################################
                      #Aligning UCE Loci
################################################################################

# make sure we are in the correct directory
cd /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/taxon-sets/all

# align the data and use edge trimming
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 198 \
    --aligner mafft \
    --cores 36 \
    --incomplete-matrix \
    --log-path log

#summary stats
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-edge-trimmed \
    --cores 36 \
    --log-path log \
    --input-format nexus \
    --output-stats unphased_align_stats.out

# run gblocks with internal trimming on the alignments
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft-nexus-edge-trimmed \
    --input-format nexus \
    --output mafft-nexus-edge-internal-trimmed-gblocks \
    --output-format nexus\
    --cores 36 \
    --log log

#summary stats for the gblock-trimmed alignments
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-edge-internal-trimmed-gblocks \
    --cores 36 \
    --log-path log \
    --output glbocks-trim-stats.out

################################################################################
                          #Alignment Cleaning
################################################################################

# make sure we are in the correct directory
cd /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/taxon-sets/all

# align the data - turn off trimming and output FASTA
phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments  mafft-nexus-edge-internal-trimmed-gblocks \
    --output mafft-nexus-edge-internal-trimmed-gblocks-clean \
    --cores 36 \
    --log-path log

################################################################################
                #Final Data matrices with % of missing data
################################################################################

# make sure we are in the correct directory
cd /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/taxon-sets/all

#50% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean \
    --taxa 108 \
    --percent 0.50 \
    --output mafft-nexus-edge-internal-trimmed-gblocks-clean-50p \
    --cores 36 \
    --log-path log

#75% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean \
    --taxa 108 \
    --percent 0.75 \
    --output mafft-nexus-edge-internal-trimmed-gblocks-clean-75p \
    --cores 36 \
    --log-path log

#90% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean \
    --taxa 108 \
    --percent 0.90 \
    --output mafft-nexus-edge-internal-trimmed-gblocks-clean-90p \
    --cores 36 \
    --log-path log

#100% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean \
    --taxa 108 \
    --percent 1.0 \
    --output mafft-nexus-edge-internal-trimmed-gblocks-clean-100p \
    --cores 36 \
    --log-path log


################################################################################
#Files for RAxML and svdquartets
############################################################

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
     --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean-75p\
     --output mafft-nexus-post-edge-internal-trimmed-gblocks-clean-75p-raxml \
     --charsets \
     --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
     --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean-75p \
     --output mafft-nexus-post-edge-internal-trimmed-gblocks-clean-75p-svdq \
     --charsets \
     --nexus \
     --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
     --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean-50p \
     --output mafft-nexus-post-edge-internal-trimmed-gblocks-clean-50p-raxml \
     --charsets \
     --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
     --alignments mafft-nexus-edge-internal-trimmed-gblocks-clean-50p \
     --output mafft-nexus-post-edge-internal-trimmed-gblocks-clean-50p-svdq \
     --charsets \
     --nexus \
     --log-path log
