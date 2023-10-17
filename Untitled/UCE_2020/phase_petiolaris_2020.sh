#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=144:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node
#SBATCH --mem=240G   # maximum memory per node
#SBATCH --job-name="phase_pet"
#SBATCH --mail-user=kevinq@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluceBio

#current working directory
cd /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/




################################################################################
                          #Extracting UCE Loci
################################################################################
# create an output directory for this taxon set - this just keeps
# things neat

#mkdir taxon-sets/petiolaris2

# create the data matrix configuration file
#phyluce_assembly_get_match_counts \
#    --locus-db uce-search-results/probe.matches.sqlite \
#    --taxon-list-config petiolaris-set.conf \
#    --taxon-group 'all' \
#    --incomplete-matrix \
#    --output taxon-sets/petiolaris2/petiolaris-taxa-incomplete.conf

# change to the taxon-sets/all directory
cd taxon-sets/petiolaris2

# make a log directory to hold our log files - this keeps things neat
#mkdir log


# get FASTA data for taxa in our taxon set
#phyluce_assembly_get_fastas_from_match_counts \
#    --contigs /work/LAS/phylo-lab/kevinq/Fig_wasp_phylogeography/Data/UCE/spades_assemblies/contigs \
#    --locus-db ../../uce-search-results/probe.matches.sqlite \
#    --match-count-output petiolaris-taxa-incomplete.conf \
#    --output petiolaris-taxa-incomplete.fasta \
#    --incomplete-matrix petiolaris-taxa-incomplete.incomplete \
#    --log-path log

# align the data but do NOT trim
#phyluce_align_seqcap_align \
#    --fasta petiolaris-taxa-incomplete.fasta \
#    --output mafft-fasta-edge-trimmed \
#    --taxa 108 \
#    --aligner mafft \
#    --cores 36 \
#    --incomplete-matrix \
#    --output-format fasta \
#    --log-path log

#explode the untrimmed alignments
#phyluce_align_explode_alignments \
#    --alignments mafft-fasta-edge-trimmed \
#    --input-format fasta \
#    --output mafft-fasta-edge-trimmed-exploded \
#    --by-taxon

#get summary stats on the FASTAS
#for i in mafft-fasta-edge-trimmed-exploded/*.fasta;
#do
#    phyluce_assembly_get_fasta_lengths --input $i --csv >> exploded.txt;
#done


#map reads to contigs
# phyluce_snp_bwa_multiple_align \
#     --config ../../phasing_petiolaris.conf \
#     --output multialign-bams-untrimmed \
#     --cores 36 \
#     --log-path log \
#     --mem

# #phase alleles
# phyluce_snp_phase_uces \
#     --config ../../phasing_petiolaris.conf \
#     --bams multialign-bams-untrimmed \
#     --cores 36 \
#     --output multialign-bams-untrimmed-phased-reads

#realign phased loci and edge trim
phyluce_align_seqcap_align \
    --fasta multialign-bams-untrimmed-phased-reads/fastas/joined_allele_sequences_all_samples.fasta \
    --output phased-mafft-fasta-post-edge-trimmed \
    --taxa 108 \
    --aligner mafft \
    --ambiguous \
    --cores 4 \
    --incomplete-matrix \
    --output-format fasta \
    --log-path log

#summary stats for the phased alignments
phyluce_align_get_align_summary_data \
    --alignments phased-mafft-fasta-post-edge-trimmed \
    --cores 36 \
    --input-format fasta \
    --log-path log

# run gblocks with internal trimming on the phased alignments
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments phased-mafft-fasta-post-edge-trimmed \
    --input-format fasta \
    --output phased-mafft-nexus-post-edge-internal-trimmed-gblocks \
    --cores 36 \
    --log log

#summary stats for the phased gblock-trimmed alignments
phyluce_align_get_align_summary_data \
    --alignments phased-mafft-nexus-post-edge-internal-trimmed-gblocks \
    --cores 36 \
    --log-path log

# align the data - turn off trimming and output FASTA
phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments phased-mafft-nexus-post-edge-internal-trimmed-gblocks \
    --input-format nexus \
    --output phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean \
    --output-format nexus \
    --cores 36 \
    --log-path log

################################################################################
                #Final Data matrices with % of missing data
################################################################################

# the integer following --taxa is the number of TOTAL taxa
# and I use "100p" to denote the 100% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean \
    --taxa 216 \
    --percent 1.0 \
    --output phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean-100p \
    --cores 36 \
    --log-path log

# the integer following --taxa is the number of TOTAL taxa
# and I use "90p" to denote the 90% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean \
    --taxa 216 \
    --percent 0.90 \
    --output phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean-90p \
    --cores 36 \
    --log-path log

# the integer following --taxa is the number of TOTAL taxa
# and I use "70p" to denote the 70% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean \
    --taxa 216 \
    --percent 0.70 \
    --output phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean-70p \
    --cores 36 \
    --log-path log

# 50% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean \
    --taxa 216 \
    --percent 0.50 \
    --output phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean-50p \
    --cores 36 \
    --log-path log

# #################################################
# ### Build input files for RAxML/svdq analysis ###
# #################################################
#
# # build the concatenated data matrix for raxml
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-100p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-100p-raxml \
#     --charsets \
#     --log-path log
#
# # build the concatenated data matrix for svdquartets
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-100p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-100p-svdq \
#     --charsets \
#     --nexus \
#     --log-path log
#
# # build the concatenated data matrix for raxml
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-90p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-90p-raxml \
#     --charsets \
#     --log-path log
#
# # build the concatenated data matrix for svdquartets
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-90p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-90p-svdq \
#     --charsets \
#     --nexus \
#     --log-path log
#
# # build the concatenated data matrix for raxml
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-70p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-70p-raxml \
#     --charsets \
#     --log-path log
#
# # build the concatenated data matrix for svdquartets
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-70p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-70p-svdq \
#     --charsets \
#     --nexus \
#     --log-path log
#
# # build the concatenated data matrix for raxml
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-50p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-50p-raxml \
#     --charsets \
#     --log-path log
#
# # build the concatenated data matrix for svdquartets
# phyluce_align_format_nexus_files_for_raxml \
#     --alignments phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-50p \
#     --output phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-50p-svdq \
#     --charsets \
#     --nexus \
#     --log-path log
